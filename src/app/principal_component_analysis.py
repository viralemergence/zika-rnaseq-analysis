from argparse import ArgumentParser
from contextlib import contextmanager, redirect_stdout, redirect_stderr
import matplotlib.pyplot as plt # type: ignore
from os import devnull, environ
import pandas as pd # type: ignore
from pydeseq2.dds import DeseqDataSet # type: ignore

environ["NUMBA_CACHE_DIR"] = "/tmp/" # Needed for scanpy to import properly
import scanpy as sc # type: ignore
from warnings import catch_warnings, simplefilter

@contextmanager
def suppress_stdout_stderr():
    with open(devnull, "w") as null_handle:
        with redirect_stdout(null_handle) as out, redirect_stderr(null_handle) as err:
            yield(out, err)

class PrincipalComponentAnalysisManager:
    def __init__(self):
        self.pca_figure_dir = "/src/data/pydeseq2/pca/"
        sc.settings.figdir = self.pca_figure_dir

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

                print("Starting pyDEseq2 analysis") # TODO: Could put all this in a "differential expression" function for easier reading
                with catch_warnings():
                    simplefilter("ignore")
                    dds = DeseqDataSet(counts=gene_counts_by_virus,
                                    metadata=sample_metadata_by_virus,
                                    design_factors=design_factors)
                with suppress_stdout_stderr(): # NOTE: Primarily suppresses redundant messages in prod, but could suppress actual errors important for dev work
                    dds.deseq2()
                
                dds.obs["Lib. Prep Batch"] = dds.obs["Lib. Prep Batch"].astype(int).astype(str)
                self.perform_principal_component_analysis(dds, cell_line, virus, pca_color_factors)

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
        plt.rcParams["svg.fonttype"] = "none"
        sc.tl.pca(dds)
        parameters = {"size": 200, "annotate_var_explained": True}
        for design_factor in design_factors:
            sc.pl.pca(dds, color=design_factor, save=f"_{cell_line}_{virus}_{design_factor}.svg", **parameters)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-c", "--counts", type=str, required=True)
    parser.add_argument("-m", "--metadata", type=str, required=True)
    args = parser.parse_args()

    gene_counts = pd.read_csv(args.counts, index_col=0)
    sample_metadata = pd.read_csv(args.metadata, index_col=0)
    
    pcam = PrincipalComponentAnalysisManager()
    pcam.run(gene_counts, sample_metadata)