from argparse import ArgumentParser
from contextlib import contextmanager, redirect_stdout, redirect_stderr
import matplotlib.pyplot as plt # type: ignore
from os import devnull, environ
import pandas as pd # type: ignore
from pathlib import Path
from pydeseq2.dds import DeseqDataSet # type: ignore

environ["NUMBA_CACHE_DIR"] = "/tmp/" # Needed for scanpy to import properly
import scanpy as sc # type: ignore
import scipy.cluster.hierarchy as sch # type: ignore
import subprocess
from typing import Union
from warnings import catch_warnings, simplefilter

@contextmanager
def suppress_stdout_stderr():
    with open(devnull, "w") as null_handle:
        with redirect_stdout(null_handle) as out, redirect_stderr(null_handle) as err:
            yield(out, err)

class DifferentialExpressionAnalysis:
    def __init__(self):
        self.pca_figure_dir = "/src/data/pydeseq2/pca/"
        sc.settings.figdir = self.pca_figure_dir
        self.dendrogram_figure_dir = "/src/data/pydeseq2/dendrogram/"
        self.lrt_dir = Path("/src/data/pydeseq2/lrt/")
        self.degpatterns_dir = Path("/src/data/pydeseq2/degpatterns/")

    def run(self, gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame) -> None:
        design_factors = ["Time", "Virus"] # Doesn't like Virus AND Treatment (redundancy probably)
        viruses = ["MR", "PRV"] # TODO: Probably better to extract from metadata file
        pca_color_factors = ["Lib. Prep Batch", "Time", "Virus"] # TODO: Probably better to extract, maybe from a yaml

        cell_lines = self.extract_cell_lines(sample_metadata)
        print(f"Metadata listed cell lines: {cell_lines}")

        for cell_line in cell_lines:
            if cell_line in ["HypNi"]: # NOTE: Will want to remove for future projects
                print(f"\nSkipping cell line: {cell_line}")
                continue
            print(f"\nStarting on cell line: {cell_line}\n----------")
            gene_counts_by_cell_line, sample_metadata_by_cell_line = self.filter_for_cell_line(gene_counts, sample_metadata, cell_line)
            
            for virus in viruses:
                print(f"\nStarting on virus: {virus}")
                gene_counts_by_virus, sample_metadata_by_virus = self.filter_for_virus(gene_counts_by_cell_line, sample_metadata_by_cell_line, virus)
                lrt_results = self.perform_likelihood_ratio_test(gene_counts_by_virus, sample_metadata_by_virus, self.lrt_dir, cell_line, virus)

                print("Starting pyDEseq2 analysis") # TODO: Could put all this in a "differential expression" function for easier reading
                with catch_warnings():
                    simplefilter("ignore")
                    dds = DeseqDataSet(counts=gene_counts_by_virus,
                                    metadata=sample_metadata_by_virus,
                                    design_factors=design_factors)
                with suppress_stdout_stderr(): # NOTE: Primarily suppresses redundant messages in prod, but could suppress actual errors important for dev work
                    dds.deseq2()
                
                gene_clusters = self.perform_degpatterns_clustering(lrt_results, dds, sample_metadata_by_virus, self.degpatterns_dir, cell_line, virus)
                continue

                self.lrt_sanity_check_graphs(lrt_results, dds)
                return

                dds.obs["Lib. Prep Batch"] = dds.obs["Lib. Prep Batch"].astype(int).astype(str)
                self.perform_principal_component_analysis(dds, cell_line, virus, pca_color_factors)
                
                self.perform_hierarchical_clustering(dds, cell_line, virus, self.dendrogram_figure_dir)

    @staticmethod
    def remove_discrepant_sample_ids(gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
        gene_count_sample_ids = set(gene_counts.index)
        metadata_sample_ids = set(sample_metadata.index)
        shared_ids = gene_count_sample_ids & metadata_sample_ids

        gene_counts = gene_counts[gene_counts.index.isin(shared_ids)]
        sample_metadata = sample_metadata[sample_metadata.index.isin(shared_ids)]
        return gene_counts, sample_metadata
    
    @staticmethod
    def extract_cell_lines(sample_metadata: pd.DataFrame) -> set[str]:
        cell_lines = set(sample_metadata["Cell Line"].to_list())
        cell_lines.discard("nan")
        return cell_lines
    
    @staticmethod
    def filter_for_cell_line(gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame, cell_line: str) -> tuple[pd.DataFrame, pd.DataFrame]:
        sample_metadata = sample_metadata[sample_metadata["Cell Line"] == cell_line]
        gene_counts = gene_counts[gene_counts.index.isin(sample_metadata.index)]
        return gene_counts, sample_metadata
    
    @staticmethod
    def filter_for_virus(gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame, virus: str) -> tuple[pd.DataFrame, pd.DataFrame]:
        virus_list = [virus, "No_Virus"]
        sample_metadata = sample_metadata[sample_metadata["Virus"].isin(virus_list)]
        gene_counts = gene_counts[gene_counts.index.isin(sample_metadata.index)]
        return gene_counts, sample_metadata

    @staticmethod
    def perform_principal_component_analysis(dds: DeseqDataSet, cell_line: str, virus: str, design_factors: list[str]) -> None:
        sc.tl.pca(dds)
        parameters = {"size": 200, "annotate_var_explained": True}
        for design_factor in design_factors:
            sc.pl.pca(dds, color=design_factor, save=f"_{cell_line}_{virus}_{design_factor}.png", **parameters)

    @staticmethod
    def perform_hierarchical_clustering(dds: DeseqDataSet, cell_line: str, virus: str, figure_dir: str) -> None:
        Z = sch.linkage(dds.obsm["X_pca"], method="complete", metric="correlation")
        fig, ax = plt.subplots()
        sch.dendrogram(Z, ax=ax, orientation="left", labels=dds.obs["DesiredFileName"].tolist())
        labels = ax.get_ymajorticklabels()
        for label in labels:
            virus_treatment = str(label).split("_")[1]
            timepoint = str(label).split("_")[2]
            if virus_treatment == "NA" or timepoint == "T0":
                color = "r"
            else:
                color = "b"
            label.set_color(color)
        plt.savefig(f"{figure_dir}dendrogram_{cell_line}_{virus}.png", bbox_inches="tight")
        
    @classmethod
    def perform_likelihood_ratio_test(cls, gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame, write_directory: Path, cell_line: str, virus: str) -> pd.DataFrame:
        gene_counts, sample_metadata = cls.remove_zero_time_point(gene_counts, sample_metadata)
        gene_count_path, sample_metadata_path = cls.write_input_data_for_log_ratio_test(gene_counts, sample_metadata, write_directory)
        r_lrt_results_path = write_directory / f"{cell_line}_{virus}_LRT_results.csv"

        cls.run_r_lrt_command(gene_count_path, sample_metadata_path, r_lrt_results_path)
        return pd.read_csv(r_lrt_results_path, index_col=0)

    @classmethod
    def remove_zero_time_point(cls, gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame) -> tuple[pd.DataFrame]:
        sample_metadata = sample_metadata[sample_metadata["Time"] != 0]
        return cls.remove_discrepant_sample_ids(gene_counts, sample_metadata)

    @staticmethod
    def write_input_data_for_log_ratio_test(gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame, directory: Path) -> tuple[str]:
        gene_count_outpath = directory / "gene_counts_R_input.csv"
        sample_metadata_outpath = directory / "sample_metadata_R_input.csv"
        gene_counts.to_csv(gene_count_outpath, index=True)
        sample_metadata.to_csv(sample_metadata_outpath, index=True)
        return gene_count_outpath, sample_metadata_outpath
    
    @staticmethod
    def run_r_lrt_command(gene_count_path: str, sample_metadata_path: str, r_lrt_results_path: Path) -> None:
        print("Starting LRT in R")
        r_lrt_command = ["Rscript", "/src/app/log_ratio_test.R",
                         gene_count_path, sample_metadata_path, r_lrt_results_path]

        p = subprocess.Popen(r_lrt_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        while p.poll() is None and (line := p.stdout.readline()) != "":
            pass
        p.wait()

        if p.poll() != 0:
            print(f"Exit code: {p.poll()}")
            for line in p.stderr.readlines():
                print(line.strip())
            raise Exception("R DESeq2 LRT did not complete successfully")

    @classmethod
    def lrt_sanity_check_graphs(cls, lrt_results: pd.DataFrame, dds: DeseqDataSet) -> None:
        normalized_counts = cls.extract_normalized_count_df_from_dds(dds)
        lrt_gene_ids_of_interest_parameters = [{"p_threshold": 0.05, "direction": "smaller", "head": True},
                                               {"p_threshold": 0.5, "direction": "bigger", "head": True}]

        for parameters in lrt_gene_ids_of_interest_parameters:
            gene_ids_of_interest = cls.extract_lrt_gene_ids_of_interest(lrt_results, **parameters)
            normalized_counts_of_interest = normalized_counts[gene_ids_of_interest].copy()
            normalized_counts_of_interest["Time"] = dds.obs["Time"].astype(float).astype(int)
            normalized_counts_of_interest["Virus"] = dds.obs["Virus"]
            normalized_counts_of_interest["Color"] = normalized_counts_of_interest["Virus"].apply(lambda x: "Orange" if x == "No-Virus" else "Blue")

            for gene_id in gene_ids_of_interest:
                ax = normalized_counts_of_interest.groupby(["Time", "Virus"])[gene_id].mean().unstack().plot(legend=True)
                normalized_counts_of_interest.plot(x="Time", y=gene_id, kind="scatter", ax=ax, color=normalized_counts_of_interest["Color"])

                # TODO: Still need to add cell name and virus
                figure_filename = f"{gene_id}_{parameters['p_threshold']}_{parameters['direction']}.png"
                figure_outpath = f"/src/data/pydeseq2/gene_counts_by_time/{figure_filename}"
                plt.savefig(figure_outpath, bbox_inches="tight")

    @staticmethod
    def extract_normalized_count_df_from_dds(dds: DeseqDataSet) -> pd.DataFrame:
        normalized_counts = dds.layers["normed_counts"]
        gene_ids = dds._var.index.to_list()
        sample_ids = dds.obsm["design_matrix"].index.to_list()
        return pd.DataFrame(normalized_counts, index=sample_ids, columns=gene_ids)
    
    @staticmethod
    def extract_lrt_gene_ids_of_interest(lrt_results: pd.DataFrame, p_threshold: float, direction: str, head: bool = False) -> list[str]:
        if direction == "smaller":
            significant_lrt_results = lrt_results[lrt_results["padj"] < p_threshold].sort_values("padj")
        if direction == "bigger":
            significant_lrt_results = lrt_results[lrt_results["padj"] > p_threshold].sort_values("padj")
        
        if head:
            return significant_lrt_results.head(5).index.tolist()
        return significant_lrt_results.index.tolist()

    @classmethod
    def perform_degpatterns_clustering(cls, lrt_results: pd.DataFrame, dds: DeseqDataSet, sample_metadata: pd.DataFrame, write_directory: Path, cell_line: str, virus: str) -> Union[dict[int], None]:
        gene_ids_of_interest = cls.extract_lrt_gene_ids_of_interest(lrt_results, 0.05, "smaller")
        if len(gene_ids_of_interest) == 0:
            print("No LRT significant padj: skipping degpatterns clustering")
            return None
        normalized_counts = cls.extract_normalized_count_df_from_dds(dds)
        normalized_counts, sample_metadata = cls.remove_zero_time_point(normalized_counts, sample_metadata) # NOTE: May not want this step
        normalized_counts_of_interest = normalized_counts[gene_ids_of_interest].copy()
        # NOTE: A stupid way to get the No_Virus category at the end of a sort so the labels are consistently colored by degpatterns:
        sample_metadata["Virus"] = sample_metadata["Virus"].apply(lambda x: "~No_Virus" if x == "No_Virus" else x)

        gene_count_path, sample_metadata_path = cls.write_input_data_for_degpattern_clustering(normalized_counts_of_interest, sample_metadata, write_directory)
        degpattern_results_path = write_directory / f"{cell_line}_{virus}_gene_clusters.csv"
        degpattern_figure_path = write_directory / f"{cell_line}_{virus}_gene_clusters.pdf"
        cls.run_r_degpatterns_command(gene_count_path, sample_metadata_path, degpattern_results_path, degpattern_figure_path)
        return dict(pd.read_csv(degpattern_results_path).values)

    @staticmethod
    def write_input_data_for_degpattern_clustering(gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame, directory: Path) -> tuple[str]:
        gene_count_outpath = directory / "gene_counts_degpattern_input.csv"
        sample_metadata_outpath = directory / "sample_metadata_degpattern_input.csv"
        gene_counts.to_csv(gene_count_outpath, index=True)
        sample_metadata.to_csv(sample_metadata_outpath, index=True)
        return gene_count_outpath, sample_metadata_outpath

    @staticmethod
    def run_r_degpatterns_command(gene_count_path: str, sample_metadata_path: str, degpattern_results_path: Path, degpattern_figure_path: Path) -> None:
        print("Starting degPatterns in R")
        r_degpatterns_command = ["Rscript", "/src/app/degpattern_clustering.R",
                         gene_count_path, sample_metadata_path, degpattern_results_path, degpattern_figure_path]

        p = subprocess.Popen(r_degpatterns_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        while p.poll() is None and (line := p.stdout.readline().strip()) != "":
            pass
        p.wait()

        if p.poll() != 0:
            print(f"Exit code: {p.poll()}")
            for line in p.stderr.readlines():
                print(line.strip())
            raise Exception("R degPatterns did not complete successfully")

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-c", "--counts", type=str, required=True)
    parser.add_argument("-m", "--metadata", type=str, required=True)
    args = parser.parse_args()

    gene_counts = pd.read_csv(args.counts, index_col=0)
    sample_metadata = pd.read_csv(args.metadata, index_col=0)
    
    dea = DifferentialExpressionAnalysis()
    dea.run(gene_counts, sample_metadata)