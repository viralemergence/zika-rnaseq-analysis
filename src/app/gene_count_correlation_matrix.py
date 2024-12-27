from argparse import ArgumentParser
from collections import defaultdict
from contextlib import contextmanager, redirect_stdout, redirect_stderr
from itertools import chain
import matplotlib.patches as mpatches # type: ignore
import matplotlib.pyplot as plt # type: ignore
from os import devnull
import pandas as pd # type: ignore
from pathlib import Path
from pydeseq2.dds import DeseqDataSet # type: ignore
import seaborn as sns # type: ignore
from warnings import catch_warnings, simplefilter

@contextmanager
def suppress_stdout_stderr():
    with open(devnull, "w") as null_handle:
        with redirect_stdout(null_handle) as out, redirect_stderr(null_handle) as err:
            yield(out, err)

class GeneCountCorrelations:
    def __init__(self, cell_lines: list[str], outdir: Path) -> None:
        self.cell_lines = cell_lines
        self.outdir = outdir

    def run(self, gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame) -> None:
        design_factors = ["Time", "Virus"]

        for cell_line in self.cell_lines:
            print(f"\nStarting on cell line: {cell_line}")
            gene_counts_by_cell_line, sample_metadata_by_cell_line = self.filter_for_cell_line(gene_counts, sample_metadata, cell_line)
            dds = self.pydeseq2_normalization(gene_counts_by_cell_line, sample_metadata_by_cell_line, design_factors)
            
            normalized_counts_by_cell_line = self.extract_normalized_count_df_from_dds(dds)
            normalized_counts_by_cell_line, sample_metadata_by_cell_line = self.remove_zero_time_point(normalized_counts_by_cell_line, sample_metadata_by_cell_line)

            normalized_counts = self.transform_normalized_counts(normalized_counts_by_cell_line, dds)
            per_gene_stats = normalized_counts.aggregate(["mean", "std"], axis=0).round(2)

            drop_genes = [] # NOTE: Turn into function
            for gene in per_gene_stats:
                if per_gene_stats.loc["mean", gene] == 0:
                    drop_genes.append(gene)
            normalized_counts = normalized_counts.drop(drop_genes, axis=1)
            
            custom_order = ["No Virus", "MR", "PRV"]
            normalized_counts = normalized_counts.sort_index(key=lambda x: x.map(dict(zip(custom_order, range(len(custom_order))))))
            
            method = "spearman"
            correlation_matrix = normalized_counts.T.corr(method=method)
            correlation_matrix = correlation_matrix.apply(lambda x: round(x, 4))
            max_value = correlation_matrix.max().max()
            min_value = correlation_matrix.min().min()
            print(f"{max_value}, {min_value}")
            
            fig, ax = plt.subplots()
        
            heatmap = sns.heatmap(correlation_matrix, ax=ax, vmin=0.9)
            cbar = heatmap.collections[0].colorbar
            cbar.set_label("Spearman Coefficient", labelpad=10)
            columns = [4, 8]
            for col in columns:
                heatmap.axvline(col, color="white", lw=3)
                heatmap.axhline(col, color="white", lw=3)
            heatmap.set_xlabel("")
            heatmap.set_ylabel("")
            
            heatmap.tick_params(left=False, bottom=False)

            figure_filename = f"{cell_line}_{method}_heatmap.png"
            figure_outpath = self.outdir / figure_filename
            fig.savefig(figure_outpath, bbox_inches="tight")
            plt.close(fig)
    
    @staticmethod
    def filter_for_cell_line(gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame, cell_line: str) -> tuple[pd.DataFrame, pd.DataFrame]:
        sample_metadata = sample_metadata[sample_metadata["Cell Line"] == cell_line]
        gene_counts = gene_counts[gene_counts.index.isin(sample_metadata.index)]
        return gene_counts, sample_metadata

    @staticmethod
    def pydeseq2_normalization(gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame, design_factors: list[str]) -> DeseqDataSet:
        print("Starting pyDEseq2 analysis")
        with catch_warnings():
            simplefilter("ignore")
            dds = DeseqDataSet(counts=gene_counts,
                            metadata=sample_metadata,
                            design_factors=design_factors)
        with suppress_stdout_stderr(): # NOTE: Primarily suppresses redundant messages in prod, but could suppress actual errors important for dev work
            dds.deseq2()
        return dds

    @staticmethod
    def extract_normalized_count_df_from_dds(dds: DeseqDataSet) -> pd.DataFrame:
        normalized_counts = dds.layers["normed_counts"]
        gene_ids = dds._var.index.to_list()
        sample_ids = dds.obsm["design_matrix"].index.to_list()
        return pd.DataFrame(normalized_counts, index=sample_ids, columns=gene_ids)

    @classmethod
    def remove_zero_time_point(cls, gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame) -> tuple[pd.DataFrame]:
        sample_metadata = sample_metadata[sample_metadata["Time"] != 0]
        return cls.remove_discrepant_sample_ids(gene_counts, sample_metadata)

    @staticmethod
    def remove_discrepant_sample_ids(gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
        gene_count_sample_ids = set(gene_counts.index)
        metadata_sample_ids = set(sample_metadata.index)
        shared_ids = gene_count_sample_ids & metadata_sample_ids

        gene_counts = gene_counts[gene_counts.index.isin(shared_ids)]
        sample_metadata = sample_metadata[sample_metadata.index.isin(shared_ids)]
        return gene_counts, sample_metadata

    @staticmethod
    def transform_normalized_counts(normalized_counts: pd.DataFrame, dds: DeseqDataSet) -> pd.DataFrame:
        normalized_counts = normalized_counts.copy()
        normalized_counts["Time"] = dds.obs["Time"].astype(float).astype(int)
        normalized_counts["Virus"] = dds.obs["Virus"]
        normalized_counts["Virus"] = normalized_counts["Virus"].apply(lambda x : x.replace("-", " "))
        grouped_counts = normalized_counts.groupby(["Virus", "Time"])
        return grouped_counts.mean().round(2)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-c", "--counts", type=str, required=True)
    parser.add_argument("-m", "--metadata", type=str, required=True)
    parser.add_argument("-l", "--cell_lines", nargs="+", required=True)
    parser.add_argument("-o", "--outdir", type=str, required=True)

    args = parser.parse_args()

    gene_counts = pd.read_csv(args.counts, index_col=0)
    sample_metadata = pd.read_csv(args.metadata, index_col=0)
    
    gcc = GeneCountCorrelations(args.cell_lines, Path(args.outdir))
    gcc.run(gene_counts, sample_metadata)