from argparse import ArgumentParser
from contextlib import contextmanager, redirect_stdout, redirect_stderr
from csv import reader
from math import log2
import matplotlib.pyplot as plt # type: ignore
from os import devnull
import pandas as pd # type: ignore
from pathlib import Path
from pydeseq2.dds import DeseqDataSet # type: ignore

from warnings import catch_warnings, simplefilter

@contextmanager
def suppress_stdout_stderr():
    with open(devnull, "w") as null_handle:
        with redirect_stdout(null_handle) as out, redirect_stderr(null_handle) as err:
            yield(out, err)

class GeneRelativeAbundance:
    def __init__(self, genes_path: Path):
        self.genes_of_interest = self.extract_genes_of_interest(genes_path)
        self.cell_line = genes_path.stem.split("_")[0]
        self.virus = genes_path.stem.split("_")[1]

    @staticmethod
    def extract_genes_of_interest(genes_path: Path, header: bool = True) -> list[str]:
        goi = []
        with genes_path.open() as inhandle:
            reader_iterator = reader(inhandle)
            if header:
                header_line = next(reader_iterator)
                for line in reader_iterator:
                    try:
                        if float(line[-1]) > 0.05:
                            continue
                    except ValueError:
                        continue
                    goi.append(line[0])
        return goi
            # return list(set(map(lambda x: x[0].replace(".", "-"), reader_iterator)))

    def run(self, gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame) -> None:
        print(len(self.genes_of_interest))
        design_factors = ["Time", "Virus"] # Doesn't like Virus AND Treatment (redundancy probably)
        print(f"Cell line: {self.cell_line}\nVirus: {self.virus}")

        gene_counts_by_cell_line, sample_metadata_by_cell_line = self.filter_for_cell_line(gene_counts, sample_metadata, self.cell_line)

        gene_counts_by_virus, sample_metadata_by_virus = self.filter_for_virus(gene_counts_by_cell_line, sample_metadata_by_cell_line, self.virus)
        gene_counts_by_virus, sample_metadata_by_virus = self.remove_zero_time_point(gene_counts_by_virus, sample_metadata_by_virus)

        print("Starting pyDEseq2 analysis") # TODO: Could put all this in a "differential expression" function for easier reading
        with catch_warnings():
            simplefilter("ignore")
            dds = DeseqDataSet(counts=gene_counts_by_virus,
                            metadata=sample_metadata_by_virus,
                            design_factors=design_factors)
        with suppress_stdout_stderr(): # NOTE: Primarily suppresses redundant messages in prod, but could suppress actual errors important for dev work
            dds.deseq2()

        normalized_counts = self.extract_normalized_count_df_from_dds(dds)
        genes_of_interest = self.filter_genes_of_interest(self.genes_of_interest, normalized_counts)
        
        normalized_counts_of_interest = normalized_counts[genes_of_interest].copy().applymap(lambda x: log2(x + 1)) # NOTE: May not want 0 correction
        normalized_counts_of_interest["Time"] = dds.obs["Time"].astype(float).astype(int)
        normalized_counts_of_interest["Virus"] = dds.obs["Virus"]
        
        grouped_counts = normalized_counts_of_interest.groupby(["Time", "Virus"])
        grouped_replicate_means = grouped_counts.mean().round(2)

        time_virus_stats = grouped_replicate_means.aggregate(["mean", "std"], axis=1).round(2)

        relative_abundance_zscores = grouped_replicate_means.copy()
        for time, virus in grouped_replicate_means.index:
            mean = time_virus_stats.loc[(time, virus), "mean"]
            std = time_virus_stats.loc[(time, virus), "std"]
            relative_abundance_zscores.loc[(time, virus)] = relative_abundance_zscores.loc[(time, virus)].apply(lambda x: (x-mean)/std)

        def percentile(percentile_threshold: float):
            def percentile_(series: pd.Series):
                return series.quantile(percentile_threshold)
            percentile_.__name__ = f"percentile_{percentile_threshold*100}"
            return percentile_

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

    @classmethod
    def gene_of_interest_graphs(cls, cell_line: str, virus: str, genes_of_interest: list[str], dds: DeseqDataSet) -> None:
        normalized_counts = cls.extract_normalized_count_df_from_dds(dds)
        genes_of_interest = cls.filter_genes_of_interest(genes_of_interest, normalized_counts)

        normalized_counts_of_interest = normalized_counts[genes_of_interest].copy()
        normalized_counts_of_interest["Time"] = dds.obs["Time"].astype(float).astype(int)
        normalized_counts_of_interest["Virus"] = dds.obs["Virus"]
        normalized_counts_of_interest["Virus"] = normalized_counts_of_interest["Virus"].apply(lambda x: "~No-Virus" if x == "No-Virus" else x)
        normalized_counts_of_interest = normalized_counts_of_interest.sort_values("Virus")
        normalized_counts_of_interest["Color"] = normalized_counts_of_interest["Virus"].apply(lambda x: "orange" if x == "~No-Virus" else "cornflowerblue")

        grouped_counts = normalized_counts_of_interest.groupby(["Time", "Virus"])
        for gene_id in genes_of_interest:
            ax = grouped_counts[gene_id].mean().unstack().plot(legend=True, color=["royalblue", "darkorange"])
            normalized_counts_of_interest.plot(x="Time", y=gene_id, kind="scatter", ax=ax, color=normalized_counts_of_interest["Color"], alpha=0.75)

            figure_filename = f"{cell_line}_{virus}_{gene_id}.png"
            figure_outpath = f"/src/data/pydeseq2/relative_gene_abundance/{cell_line}_{virus}/{figure_filename}"
            plt.savefig(figure_outpath, bbox_inches="tight")
            plt.close()

    @staticmethod
    def extract_normalized_count_df_from_dds(dds: DeseqDataSet) -> pd.DataFrame:
        normalized_counts = dds.layers["normed_counts"]
        gene_ids = dds._var.index.to_list()
        sample_ids = dds.obsm["design_matrix"].index.to_list()
        return pd.DataFrame(normalized_counts, index=sample_ids, columns=gene_ids)

    @staticmethod
    def filter_genes_of_interest(genes_of_interest: list[str], normalized_counts: pd.DataFrame) -> list[str]:
        filtered_genes_of_interest = []
        missing_genes_of_interest = []
        for gene in genes_of_interest:
            if gene in normalized_counts.columns:
                filtered_genes_of_interest.append(gene)
            else:
                missing_genes_of_interest.append(gene)
        print(f"Missing the following GOI: {missing_genes_of_interest}")
        return filtered_genes_of_interest

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-c", "--counts", type=str, required=True)
    parser.add_argument("-m", "--metadata", type=str, required=True)
    parser.add_argument("-g", "--genes", type=str, required=True)
    args = parser.parse_args()

    gene_counts = pd.read_csv(args.counts, index_col=0)
    sample_metadata = pd.read_csv(args.metadata, index_col=0)
    
    gra = GeneRelativeAbundance(Path(args.genes))
    gra.run(gene_counts, sample_metadata)