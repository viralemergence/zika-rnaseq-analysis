from argparse import ArgumentParser
from csv import reader
import matplotlib.pyplot as plt # type: ignore
import os
import pandas as pd # type: ignore
from pathlib import Path
import pickle
from pydeseq2.dds import DeseqDataSet # type: ignore

os.environ["NUMBA_CACHE_DIR"] = "/tmp/" # Needed for scanpy to import properly
import scanpy as sc # type: ignore
import scipy.cluster.hierarchy as sch # type: ignore
import subprocess

class SampleMetaDataManager:
    def __init__(self, metadata_path: Path) -> None:
        self.metadata_path = metadata_path
        
    def run(self) -> None:
        sample_metadata = self.extract_sample_metadata(self.metadata_path)
        sample_metadata = self.remove_rows_without_id(sample_metadata)
        sample_metadata = self.set_gene_count_id(sample_metadata)
        sample_metadata = self.set_sample_id_as_index(sample_metadata)
        sample_metadata = self.sort_by_index(sample_metadata)
        self.sample_metadata = sample_metadata
        
        self.count_id_conversion = self.set_count_id_to_sample_id_dict(self.sample_metadata)

    @staticmethod
    def extract_sample_metadata(metadata_path: Path) -> pd.DataFrame:
        return pd.read_csv(metadata_path)

    @staticmethod
    def remove_rows_without_id(sample_metadata: pd.DataFrame) -> pd.DataFrame:
        return sample_metadata.dropna(subset=["Sample ID"])
    
    @classmethod
    def set_gene_count_id(cls, sample_metadata: pd.DataFrame) -> pd.DataFrame:
        sample_metadata["gene_count_id"] = sample_metadata.apply(cls.convert_raw_id_to_count_id, axis=1)
        return sample_metadata
    
    @staticmethod
    def convert_raw_id_to_count_id(row: pd.DataFrame) -> str:
        raw_file_name = row["Raw_file_R1"]
        old_suffix = "_L001_R1_001.fastq.gz"
        count_id = raw_file_name.replace(old_suffix, "")
        return count_id
    
    @staticmethod
    def set_sample_id_as_index(sample_metadata: pd.DataFrame) -> pd.DataFrame:
        return sample_metadata.set_index("Sample ID")
    
    @staticmethod
    def sort_by_index(sample_metadata: pd.DataFrame) -> pd.DataFrame:
        return sample_metadata.sort_index()
    
    @staticmethod
    def set_count_id_to_sample_id_dict(sample_metadata: pd.DataFrame) -> dict[str]:
        gene_count_ids = sample_metadata["gene_count_id"].tolist()
        sample_ids = sample_metadata.index.tolist()
        return {count_id: sample_id for count_id, sample_id in zip(gene_count_ids, sample_ids)}

class GeneCountManager:
    def __init__(self, gene_counts_path: Path) -> None:
        self.gene_counts_path = gene_counts_path

    def run(self, count_id_conversion: dict[str], samples_to_combine_path: Path = False) -> None:
        gene_counts = self.extract_gene_counts(self.gene_counts_path)
        gene_counts = self.remove_blacklist_samples(gene_counts)
        if samples_to_combine_path:
            gene_counts = self.combine_samples(gene_counts, samples_to_combine_path)
        gene_counts = self.rename_gene_count_ids(gene_counts, count_id_conversion)
        gene_counts = self.set_gene_id_as_index(gene_counts)
        gene_counts = self.remove_low_count_rows(gene_counts)
        gene_counts = self.transpose_gene_counts(gene_counts)
        self.gene_counts = self.sort_by_index(gene_counts)

    @staticmethod
    def extract_gene_counts(gene_counts_path: Path) -> pd.DataFrame:
        return pd.read_csv(gene_counts_path)
    
    @staticmethod
    def remove_blacklist_samples(gene_counts: pd.DataFrame) -> pd.DataFrame:
        blacklist_samples = ["HypNi_ZIKV_PRVABC59_24_a_S35"]
        return gene_counts.drop(blacklist_samples, axis=1)

    @classmethod
    def combine_samples(cls, gene_counts: pd.DataFrame, samples_to_combine_path: Path) -> pd.DataFrame:
        samples_to_combine = cls.extract_samples_to_combine(samples_to_combine_path)
        for sample_1, sample_2 in samples_to_combine:
            gene_counts[sample_1] = gene_counts[sample_1] + gene_counts[sample_2]
            gene_counts = gene_counts.drop([sample_2], axis=1)
        return gene_counts

    @staticmethod
    def extract_samples_to_combine(samples_to_combine_path: Path) -> list[list]:
        with samples_to_combine_path.open() as inhandle:
            reader_iterator = reader(inhandle, delimiter=",")
            samples_to_combine = [line for line in reader_iterator]
        return samples_to_combine

    @staticmethod
    def rename_gene_count_ids(gene_counts: pd.DataFrame, id_conversion_dict: dict[str]) -> pd.DataFrame:
        return gene_counts.rename(columns=id_conversion_dict)

    @staticmethod
    def set_gene_id_as_index(gene_counts: pd.DataFrame) -> pd.DataFrame:
        gene_counts = gene_counts.rename(columns={"sample_name": "gene_id"})
        return gene_counts.set_index("gene_id")

    @staticmethod
    def remove_low_count_rows(gene_counts: pd.DataFrame, threshold: int = 10) -> pd.DataFrame:
        return gene_counts[gene_counts.sum(axis = 1) > threshold]
    
    @staticmethod
    def transpose_gene_counts(gene_counts: pd.DataFrame) -> pd.DataFrame:
        return gene_counts.T

    @staticmethod
    def sort_by_index(gene_counts: pd.DataFrame) -> pd.DataFrame:
        return gene_counts.sort_index()

class DifferentialExpressionAnalysis:
    def __init__(self):
        self.pca_figure_dir = "/src/data/pydeseq2/pca/"
        sc.settings.figdir = self.pca_figure_dir
        self.dendrogram_figure_dir = "/src/data/pydeseq2/dendrogram/"
        self.lrt_dir = "/src/data/pydeseq2/lrt/"
        self.degpatterns_dir = Path("/src/data/pydeseq2/degpatterns/")

    def run(self, gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame) -> None:
        design_factors = ["Time", "Virus"] # Doesn't like Virus AND Treatment (redundancy probably)
        viruses = ["MR", "PRV"]
        pca_color_factors = ["Lib. Prep Batch", "Time", "Virus"]

        gene_counts, sample_metadata = self.remove_discrepant_sample_ids(gene_counts, sample_metadata)
        cell_lines = self.extract_cell_lines(sample_metadata)
        print(cell_lines)

        for cell_line in cell_lines:
            if cell_line != "R06E":
                continue
            print(f"Starting on cell line: {cell_line}")
            gene_counts_by_cell_line, sample_metadata_by_cell_line = self.filter_for_cell_line(gene_counts, sample_metadata, cell_line)
            gene_counts_by_cell_line, sample_metadata_by_cell_line = self.pickle_and_rectify_batch_effect(gene_counts_by_cell_line, sample_metadata_by_cell_line, cell_line)
            
            for virus in viruses:
                gene_counts_by_virus, sample_metadata_by_virus = self.filter_for_virus(gene_counts_by_cell_line, sample_metadata_by_cell_line, virus)
                lrt_results = self.perform_likelihood_ratio_test(gene_counts_by_virus, sample_metadata_by_virus, self.lrt_dir)

                dds = DeseqDataSet(counts=gene_counts_by_virus,
                                metadata=sample_metadata_by_virus,
                                design_factors=design_factors)
                dds.deseq2()
                
                gene_clusters = self.perform_degpatterns_clustering(lrt_results, dds, sample_metadata_by_virus, self.degpatterns_dir)
                print(gene_clusters)
                return

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
    
    @classmethod
    def pickle_and_rectify_batch_effect(cls, gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame, cell_line: str) -> tuple[pd.DataFrame, pd.DataFrame]:
        gene_counts_pickled_path = Path(f"/src/data/pydeseq2/pickles/{cell_line}_gene_counts.pkl")
        sample_metadata_pickled_path = Path(f"/src/data/pydeseq2/pickles/{cell_line}_sample_metadata.pkl")
        if gene_counts_pickled_path.is_file():
            print(f"{cell_line} pickle detected. Loading now")
            with gene_counts_pickled_path.open("rb") as inhandle:
                gene_counts = pickle.load(inhandle)
            with sample_metadata_pickled_path.open("rb") as inhandle:
                sample_metadata = pickle.load(inhandle)
        else:
            gene_counts, sample_metadata = cls.rectify_batch_effect(gene_counts, sample_metadata)
            with gene_counts_pickled_path.open("wb") as outhandle:
                pickle.dump(gene_counts, outhandle)
            with sample_metadata_pickled_path.open("wb") as outhandle:
                pickle.dump(sample_metadata, outhandle)
        return gene_counts, sample_metadata
    
    @classmethod
    def rectify_batch_effect(cls, gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
        from inmoose.pycombat import pycombat_seq # type: ignore # Need to import inside function to prevent module overlap issues
        print("Starting batch effect correction (this may take awhile)")
        batches = sample_metadata["Lib. Prep Batch"].to_list()
        covariates = sample_metadata[["Time", "Virus"]]
        gene_counts = pycombat_seq(gene_counts.T, batches, covar_mod=covariates).T.astype("int")
        return cls.correct_negative_numbers(gene_counts), sample_metadata
    
    @staticmethod
    def correct_negative_numbers(gene_counts: pd.DataFrame) -> pd.DataFrame:
        negative_count = (gene_counts < 0).sum().sum()
        if negative_count > 0:
            print(f"\nWARNING: {negative_count} negative values detected after batch effect correction")
            gene_counts[gene_counts < 0] = 0
            print("Negative values have been set to 0\n")
            return gene_counts
        return gene_counts
    
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
    def perform_likelihood_ratio_test(cls, gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame, write_directory: Path) -> pd.DataFrame:
        gene_counts, sample_metadata = cls.remove_zero_time_point(gene_counts, sample_metadata)
        gene_count_path, sample_metadata_path = cls.write_input_data_for_log_ratio_test(gene_counts, sample_metadata, write_directory)
        r_lrt_results_path = Path(f"{write_directory}LRT_results.csv")

        cls.run_r_lrt_command(gene_count_path, sample_metadata_path, r_lrt_results_path)
        return pd.read_csv(r_lrt_results_path, index_col=0)

    @classmethod
    def remove_zero_time_point(cls, gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame) -> tuple[pd.DataFrame]:
        sample_metadata = sample_metadata[sample_metadata["Time"] != 0]
        return cls.remove_discrepant_sample_ids(gene_counts, sample_metadata)

    @staticmethod
    def write_input_data_for_log_ratio_test(gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame, directory: Path) -> tuple[str]:
        gene_count_outpath = f"{directory}gene_counts_R_input.csv"
        sample_metadata_outpath = f"{directory}sample_metadata_R_input.csv"
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
        print(f"Exit code: {p.poll()}")

        if p.poll() != 0:
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
    def perform_degpatterns_clustering(cls, lrt_results: pd.DataFrame, dds: DeseqDataSet, sample_metadata: pd.DataFrame, write_directory: Path) -> dict[int]:
        gene_ids_of_interest = cls.extract_lrt_gene_ids_of_interest(lrt_results, 0.05, "smaller")
        normalized_counts = cls.extract_normalized_count_df_from_dds(dds)
        normalized_counts, sample_metadata = cls.remove_zero_time_point(normalized_counts, sample_metadata) # NOTE: May not want this step
        normalized_counts_of_interest = normalized_counts[gene_ids_of_interest].copy()

        gene_count_path, sample_metadata_path = cls.write_input_data_for_degpattern_clustering(normalized_counts_of_interest, sample_metadata, write_directory)
        degpattern_results_path = write_directory / "gene_clusters.csv"
        cls.run_r_degpatterns_command(gene_count_path, sample_metadata_path, degpattern_results_path)
        return dict(pd.read_csv(degpattern_results_path).values)

    @staticmethod
    def write_input_data_for_degpattern_clustering(gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame, directory: Path) -> tuple[str]:
        gene_count_outpath = directory / "gene_counts_degpattern_input.csv"
        sample_metadata_outpath = directory / "sample_metadata_degpattern_input.csv"
        gene_counts.to_csv(gene_count_outpath, index=True)
        sample_metadata.to_csv(sample_metadata_outpath, index=True)
        return gene_count_outpath, sample_metadata_outpath

    @staticmethod
    def run_r_degpatterns_command(gene_count_path: str, sample_metadata_path: str, degpattern_results_path: Path) -> None:
        print("Starting degPatterns in R")
        r_degpatterns_command = ["Rscript", "/src/app/degpattern_clustering.R",
                         gene_count_path, sample_metadata_path, degpattern_results_path]

        p = subprocess.Popen(r_degpatterns_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        while p.poll() is None and (line := p.stdout.readline()) != "":
            pass
        p.wait()
        print(f"Exit code: {p.poll()}")

        if p.poll() != 0:
            for line in p.stderr.readlines():
                print(line.strip())
            raise Exception("R degPatterns did not complete successfully")

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-c", "--counts", type=str, required=True)
    parser.add_argument("-m", "--metadata", type=str, required=True)
    parser.add_argument("-k", "--combine", type=str, required=False)
    args = parser.parse_args()

    smdm = SampleMetaDataManager(Path(args.metadata))
    smdm.run()

    gcm = GeneCountManager(Path(args.counts))
    gcm_parameters = {"count_id_conversion": smdm.count_id_conversion}
    if args.combine:
        gcm_parameters.update({"samples_to_combine_path": Path(args.combine)})
    gcm.run(**gcm_parameters)
    
    dea = DifferentialExpressionAnalysis()
    dea.run(gcm.gene_counts, smdm.sample_metadata)