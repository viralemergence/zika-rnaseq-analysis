from argparse import ArgumentParser
from collections import defaultdict
from contextlib import contextmanager, redirect_stdout, redirect_stderr
from csv import reader
from math import ceil, log2
import matplotlib.patches as mpatches # type: ignore
import matplotlib.pyplot as plt # type: ignore
from numpy import random # type: ignore
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
    def __init__(self, genes_path: Path, cell_line: str, virus_contrast: str, outdir: Path) -> None:
        self.genes_of_interest = self.extract_genes_of_interest(genes_path)
        self.cell_line = cell_line
        self.virus_contrast = virus_contrast
        self.treatment, self.control = self.extract_conditions_from_contrast(virus_contrast)
        self.outdir = outdir

    @staticmethod
    def extract_genes_of_interest(genes_path: Path, header: bool = True) -> dict[set[str]]:
        groups = defaultdict(set)
        with genes_path.open() as inhandle:
            reader_iterator = reader(inhandle)
            if header:
                header_line = next(reader_iterator)
            for line in reader_iterator:
                groups[str(line[1])].add(line[0])
        return groups
    
    @staticmethod
    def extract_conditions_from_contrast(contrast: str) -> tuple[str, str]:
        conditions = contrast.split("_vs_")
        return conditions[0], conditions[1]

    def run(self, gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame) -> None:
        design_factors = ["Time", "Virus"]
        print(f"Cell line: {self.cell_line}\nContrast: {self.virus_contrast}")

        gene_counts_by_cell_line, sample_metadata_by_cell_line = self.filter_for_cell_line(gene_counts, sample_metadata, self.cell_line)
        gene_counts_by_virus, sample_metadata_by_virus = self.filter_for_virus_contrasts(gene_counts_by_cell_line, sample_metadata_by_cell_line,
                                                                                         self.treatment, self.control)
        
        dds = self.pydeseq2_normalization(gene_counts_by_virus, sample_metadata_by_virus, design_factors)

        normalized_counts_by_virus = self.extract_normalized_count_df_from_dds(dds)
        normalized_counts_by_virus, sample_metadata_by_virus = self.remove_zero_time_point(normalized_counts_by_virus, sample_metadata_by_virus)
        
        desired_subplots = len(self.genes_of_interest)
        subplot_columns = 3
        subplot_rows = int(ceil(desired_subplots / subplot_columns))
        combined_fig, combined_ax = plt.subplots(subplot_rows, subplot_columns, figsize=(6.4, 1.2)) # NOTE: Will want to change figsize or make arg
        #combined_fig, combined_ax = plt.subplots(subplot_rows, subplot_columns, figsize=(6.4, 2.4)) # NOTE: Will want to change figsize or make arg
        #combined_fig, combined_ax = plt.subplots(subplot_rows, subplot_columns, figsize=(6.4, 4.8)) # NOTE: Will want to change figsize or make arg
        #combined_fig, combined_ax = plt.subplots(subplot_rows, subplot_columns, figsize=(6.4, 7.2)) # NOTE: Will want to change figsize or make arg
        #combined_fig, combined_ax = plt.subplots(subplot_rows, subplot_columns, figsize=(6.4, 9.6)) # NOTE: Will want to change figsize or make arg
        for i, (group, genes) in enumerate(self.genes_of_interest.items()):
            subplot_row = i // subplot_columns
            subplot_column = i % subplot_columns

            print(f"\nStarting on group: {group}")
            genes_of_interest = self.filter_genes_of_interest(genes, normalized_counts_by_virus)

            normalized_counts_of_interest = self.extract_transform_normalized_counts_of_interest(normalized_counts_by_virus, genes_of_interest, dds)            
            per_gene_stats = normalized_counts_of_interest.aggregate(["mean", "std"], axis=0).round(2)
            
            drop_genes = [] # NOTE: Turn into function
            for gene in per_gene_stats:
                if per_gene_stats.loc["mean", gene] == 0:
                    drop_genes.append(gene)
            normalized_counts_of_interest = normalized_counts_of_interest.drop(drop_genes, axis=1)

            gene_relative_abundance_zscores = self.calculate_gene_relative_abundance_zscores(normalized_counts_of_interest, per_gene_stats)
            zscore_stats = gene_relative_abundance_zscores.aggregate(["median", "std", "min", "max",
                                                                      self.percentile(0.25), self.percentile(0.50), self.percentile(0.75)], axis=1).round(2)

            self.graph_gene_relative_abundance(gene_relative_abundance_zscores, zscore_stats, group, genes, self.cell_line, self.virus_contrast, self.outdir)
            self.combined_graph_gene_relative_abundance(gene_relative_abundance_zscores, zscore_stats, group, genes, combined_fig, combined_ax,
                                                        subplot_row, subplot_column, subplot_rows)
            
        total_subplots = subplot_columns * subplot_rows
        i += 1
        while i < total_subplots:
            subplot_row = i // subplot_columns
            subplot_column = i % subplot_columns
            combined_ax[subplot_row, subplot_column].remove()
            i += 1
        
        figure_filename = f"{self.cell_line}_{self.virus_contrast}_combined.png"
        figure_outpath = self.outdir / figure_filename
        combined_fig.savefig(figure_outpath, bbox_inches="tight")
        plt.close(combined_fig)
    
    @staticmethod
    def filter_for_cell_line(gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame, cell_line: str) -> tuple[pd.DataFrame, pd.DataFrame]:
        sample_metadata = sample_metadata[sample_metadata["Cell Line"] == cell_line]
        gene_counts = gene_counts[gene_counts.index.isin(sample_metadata.index)]
        return gene_counts, sample_metadata
    
    @staticmethod
    def filter_for_virus_contrasts(gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame, treatment: str, control: str) -> tuple[pd.DataFrame, pd.DataFrame]:
        virus_list = [treatment, control]
        sample_metadata = sample_metadata[sample_metadata["Virus"].isin(virus_list)]
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

    @staticmethod
    def extract_transform_normalized_counts_of_interest(normalized_counts: pd.DataFrame, genes_of_interest: list[str], dds: DeseqDataSet) -> pd.DataFrame:
        normalized_counts_of_interest = normalized_counts[genes_of_interest].copy().applymap(lambda x: log2(x + 1)) # NOTE: May not want 0 correction
        normalized_counts_of_interest["Time"] = dds.obs["Time"].astype(float).astype(int)
        normalized_counts_of_interest["Virus"] = dds.obs["Virus"]
        normalized_counts_of_interest["Virus"] = normalized_counts_of_interest["Virus"].apply(lambda x: "~No-Virus" if x == "No-Virus" else x)
        grouped_counts = normalized_counts_of_interest.groupby(["Time", "Virus"])
        return grouped_counts.mean().round(2)

    @staticmethod
    def calculate_gene_relative_abundance_zscores(normalized_counts_of_interest: pd.DataFrame, per_gene_stats: pd.DataFrame):
        gene_relative_abundance_zscores = normalized_counts_of_interest.copy()
        for gene in normalized_counts_of_interest.keys():
            mean = per_gene_stats.loc["mean", gene]
            std = per_gene_stats.loc["std", gene]
            gene_relative_abundance_zscores[gene] = gene_relative_abundance_zscores[gene].apply(lambda x: round((x-mean)/std, 2))
        return gene_relative_abundance_zscores

    @staticmethod
    def percentile(percentile_threshold: float):
        def percentile_(series: pd.Series) -> float:
            return series.quantile(percentile_threshold)
        percentile_.__name__ = f"percentile_{percentile_threshold*100}"
        return percentile_

    @staticmethod
    def graph_gene_relative_abundance(gene_relative_abundance_zscores: pd.DataFrame, zscore_stats: pd.DataFrame,
                                      group: str, genes: list[str], cell_line: str, virus_contrast: str, outdir: Path) -> None:
        boxplot_data = defaultdict(lambda: defaultdict(list))
        times = set()
        for time, virus in zscore_stats.index:
            quartiles = [
                zscore_stats.loc[(time, virus), "min"],
                zscore_stats.loc[(time, virus), "percentile_25.0"],
                zscore_stats.loc[(time, virus), "percentile_50.0"],
                zscore_stats.loc[(time, virus), "percentile_75.0"],
                zscore_stats.loc[(time, virus), "max"]
                ]
            boxplot_data[virus]["quartiles"].append(quartiles)
            boxplot_data[virus]["line"].append(zscore_stats.loc[(time, virus), "percentile_50.0"])
            times.add(time)
            boxplot_data[virus]["z_scores"].append(list(gene_relative_abundance_zscores.loc[(time, virus),]))
            
        fig, ax = plt.subplots()
        base_colors = ["red", "blue"]
        colors = iter(base_colors)
        for virus, data in boxplot_data.items():
            color = next(colors)
            ax.boxplot(data["quartiles"],
                        patch_artist=True,
                        boxprops=dict(facecolor="none", color=color),
                        medianprops=dict(color=color),
                        whiskerprops=dict(color=color),
                        capprops=dict(color="none"),
                        flierprops=dict(color="none", markeredgecolor="none")
                        )
            ax.plot(list(range(1, len(data["line"])+1)), data["line"], color=color)
            
            for i, z_scores in enumerate(data["z_scores"], start=1):
                x_values = random.normal(i, 0.075, size=len(z_scores))
                ax.plot(x_values, z_scores, color=color, linestyle="None", marker="o", alpha=0.2)

        ax.set_xticks(ticks=list(range(1, len(times)+1)), labels=sorted(times))
        ax.set_xlabel("Time (hr)")
        ax.set_ylabel("Z score")
        ax.set_title(f"Group: {group}, Genes: {len(genes)}")
        ax.set_ylim(bottom=-2.5, top=2.5)
        
        legend_handles = [mpatches.Patch(color=color, label=virus.replace("~","")) for color, virus in zip(base_colors, boxplot_data)]
        ax.legend(handles=legend_handles)

        ax.grid()
        ax.tick_params(grid_color="grey", grid_alpha=0.2)

        figure_filename = f"{cell_line}_{virus_contrast}_group_{group}.png"
        figure_outpath = outdir / figure_filename
        fig.savefig(figure_outpath, bbox_inches="tight")
        plt.close(fig)

    @staticmethod
    def combined_graph_gene_relative_abundance(gene_relative_abundance_zscores: pd.DataFrame, zscore_stats: pd.DataFrame,
                                      group: str, genes: list[str], fig, ax, subplot_row: int, subplot_column: int, max_rows: int) -> None:
        if max_rows == 1:
            subplot_ax = ax[subplot_column]
        else:
            subplot_ax = ax[subplot_row, subplot_column]
        
        boxplot_data = defaultdict(lambda: defaultdict(list))
        times = set()
        for time, virus in zscore_stats.index:
            quartiles = [
                zscore_stats.loc[(time, virus), "min"],
                zscore_stats.loc[(time, virus), "percentile_25.0"],
                zscore_stats.loc[(time, virus), "percentile_50.0"],
                zscore_stats.loc[(time, virus), "percentile_75.0"],
                zscore_stats.loc[(time, virus), "max"]
                ]
            boxplot_data[virus]["quartiles"].append(quartiles)
            boxplot_data[virus]["line"].append(zscore_stats.loc[(time, virus), "percentile_50.0"])
            times.add(time)
            boxplot_data[virus]["z_scores"].append(list(gene_relative_abundance_zscores.loc[(time, virus),]))
            
        base_colors = ["red", "blue"]
        colors = iter(base_colors)
        for virus, data in boxplot_data.items():
            color = next(colors)
            subplot_ax.boxplot(data["quartiles"],
                        patch_artist=True,
                        boxprops=dict(facecolor="none", color=color),
                        medianprops=dict(color=color),
                        whiskerprops=dict(color=color),
                        capprops=dict(color="none"),
                        flierprops=dict(color="none", markeredgecolor="none")
                        )
            subplot_ax.plot(list(range(1, len(data["line"])+1)), data["line"], color=color)
            
            for i, z_scores in enumerate(data["z_scores"], start=1):
                x_values = random.normal(i, 0.075, size=len(z_scores))
                subplot_ax.plot(x_values, z_scores, color=color, linestyle="None", marker="o", markersize=3, alpha=0.2)

        subplot_ax.set_xticks(ticks=list(range(1, len(times)+1)), labels=sorted(times))
        subplot_ax.set_ylim(bottom=-2.5, top=2.5)

        subplot_ax.set_title(f"Group: {group}, Genes: {len(genes)}", loc="center", fontsize=10, pad=6,
                                                          bbox={"facecolor": "lightgrey", "boxstyle":"square,pad=0.3"})
        
        if subplot_column != 0:
            subplot_ax.tick_params(left=False, labelleft=False)
        else:
            subplot_ax.set_ylabel("Z score")
            
        if (subplot_row + 1) != max_rows:
            subplot_ax.tick_params(bottom=False, labelbottom=False)
        else:
            subplot_ax.set_xlabel("Time (hr)")
            
        if subplot_column != 0 and (subplot_row + 1) == max_rows:
            legend_handles = [mpatches.Patch(color=color, label=virus.replace("~","")) for color, virus in zip(base_colors, boxplot_data)]
            fig.legend(handles=legend_handles, bbox_to_anchor=(0.9,0.5), loc="center left")
            
        subplot_ax.grid()
        subplot_ax.tick_params(grid_color="grey", grid_alpha=0.2)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-c", "--counts", type=str, required=True)
    parser.add_argument("-m", "--metadata", type=str, required=True)
    parser.add_argument("-g", "--genes", type=str, required=True)
    parser.add_argument("-l", "--cell_line", type=str, required=True)
    parser.add_argument("-v", "--virus_contrast", type=str, required=True)
    parser.add_argument("-o", "--outdir", type=str, required=True)

    args = parser.parse_args()

    gene_counts = pd.read_csv(args.counts, index_col=0)
    sample_metadata = pd.read_csv(args.metadata, index_col=0)
    
    gra = GeneRelativeAbundance(Path(args.genes), args.cell_line, args.virus_contrast, Path(args.outdir))
    gra.run(gene_counts, sample_metadata)