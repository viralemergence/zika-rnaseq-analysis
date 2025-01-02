from argparse import ArgumentParser
from collections import defaultdict
from contextlib import contextmanager, redirect_stdout, redirect_stderr
from gseapy import prerank # type: ignore
import matplotlib.pyplot as plt # type: ignore
from numpy import log10 # type: ignore
from os import devnull
import pandas as pd # type: ignore
from pathlib import Path
import seaborn as sns # type: ignore

@contextmanager
def suppress_stdout_stderr():
    with open(devnull, "w") as null_handle:
        with redirect_stdout(null_handle) as out, redirect_stderr(null_handle) as err:
            yield(out, err)

class GeneSetEnrichmentAnalysis:
    def __init__(self, count_contrasts: Path, outdir: Path) -> None:
        self.count_contrasts = self.extract_count_contrasts(count_contrasts)
        self.outdir = outdir

    @staticmethod
    def extract_count_contrasts(count_contrasts: Path) -> pd.DataFrame:
        df = pd.read_csv(count_contrasts, index_col=0)
        columns = [eval(column) for column in df.columns]
        df.columns = pd.MultiIndex.from_tuples(columns)
        return df

    def run(self) -> None:
        time_points_of_interest = [6.0, 12.0, 24.0, 48.0]
        expression_directionality = ["up", "down"]

        for contrast, deseq_results in self.count_contrasts.groupby(level=0, axis=1):
            print(f"\nStarting on contrast: {contrast}")
            deseq_results = deseq_results.droplevel(0, axis=1)
            deseq_results.columns = deseq_results.columns.set_levels(deseq_results.columns.levels[0].astype(float), level=0)
            
            gsea_collated_results = defaultdict(lambda: dict())
            
            for direction in expression_directionality:
                print(f"Starting on {direction} regulated gene sets")
                for time_point in time_points_of_interest:
                    print(f"Starting on time: {time_point}")
                    time_point_results = deseq_results[time_point].copy()
                    ranked_genes = self.rank_genes(time_point_results, direction)
                    self.run_gsea_and_append(ranked_genes, gsea_collated_results, direction, time_point)
                
            print("Starting top gene set parsing")
            top_gene_sets = self.extract_top_gene_sets(expression_directionality, time_points_of_interest, gsea_collated_results)

            print("Starting top gene set concat")
            heatmap_dataframe = self.generate_heatmap_dataframe(gsea_collated_results, time_points_of_interest, top_gene_sets)
            
            self.generate_heatmap(heatmap_dataframe, contrast, self.outdir)

    @staticmethod
    def rank_genes(time_point_results: pd.DataFrame, direction: str) -> pd.DataFrame:
        with suppress_stdout_stderr():
            if direction == "up":
                time_point_results["Rank"] = -log10(time_point_results["padj"]*time_point_results["log2FoldChange"])
            if direction == "down":
                time_point_results["Rank"] = -log10(time_point_results["padj"]*(-1*time_point_results["log2FoldChange"]))
        time_point_results = time_point_results.sort_values("Rank", ascending = False)
        time_point_results["Gene"] = time_point_results.index
        return time_point_results[["Gene", "Rank"]].reset_index(drop=True)

    @staticmethod
    def run_gsea_and_append(ranked_genes: pd.DataFrame, gsea_collated_results: dict[dict[pd.DataFrame]], direction: str, time_point: float) -> None:
        with suppress_stdout_stderr():
            gsea_results = prerank(rnk=ranked_genes, gene_sets="KEGG_2021_Human", seed=7).results
        parsed_results = []
        for term, data in gsea_results.items():
            parsed_results.append([term,
                        data["fdr"],
                        data["es"],
                        data["nes"]])
        parsed_results = pd.DataFrame(parsed_results, columns=["Term", "FDR", "ES", "NES"]).sort_values("FDR").reset_index(drop=True)
        parsed_results = parsed_results.set_index("Term")

        gsea_collated_results[direction][time_point] = parsed_results

    @staticmethod
    def extract_top_gene_sets(expression_directionality: list[str], time_points_of_interest: list[float],
                              gsea_collated_results: dict[dict[pd.DataFrame]]) -> dict[list[str]]:
        top_gene_sets = {direction: list() for direction in expression_directionality}
        first_last_time_points = [time_points_of_interest[0], time_points_of_interest[-1]]
        for direction, data in gsea_collated_results.items():
            for time_point in first_last_time_points:
                top_gene_sets_subset = list(data[time_point].head(5).index)
                for gene_set in top_gene_sets_subset:
                    if gene_set not in top_gene_sets[direction]:
                        top_gene_sets[direction].append(gene_set)
        return top_gene_sets

    @staticmethod
    def generate_heatmap_dataframe(gsea_collated_results: dict[dict[pd.DataFrame]],
                                   time_points_of_interest: list[float],
                                   top_gene_sets: dict[list[str]]) -> pd.DataFrame:
        heatmap_data = []
        for direction, gsea_data in gsea_collated_results.items():
            for time_point in time_points_of_interest:
                df = gsea_data[time_point].reindex(top_gene_sets[direction])
                df = pd.DataFrame(df.loc[top_gene_sets[direction], "NES"])
                if direction == "down":
                    df["NES"] = -1*df["NES"]
                df = df.rename(columns={"NES": f"{time_point}"})
                heatmap_data.append(df)
        heatmap_dataframe = pd.concat(heatmap_data, axis=1)
        heatmap_dataframe.columns = heatmap_dataframe.columns.astype(float)
        return heatmap_dataframe.groupby(level=0, axis=1).max()

    @staticmethod
    def generate_heatmap(heatmap_dataframe: pd.DataFrame, contrast: str, outdir: Path) -> None:
        fig, ax = plt.subplots()
    
        heatmap = sns.heatmap(heatmap_dataframe, ax=ax, vmin=-3, vmax=3, cmap="RdBu_r", square=True)

        # Title
        title = contrast.split("_vs_")
        title.insert(1, "vs")
        title = "\n".join(title)
        heatmap.set_title(title)

        # Colorbar and null values
        cbar = heatmap.collections[0].colorbar
        cbar.set_label("Normalized Enrichment Score", labelpad=10)
        heatmap.collections[0].cmap.set_bad("grey")

        # Truncating and setting labels
        heatmap.set_xlabel("")
        heatmap.set_ylabel("")
        
        labels = []
        for label in heatmap.get_yticklabels():
            text = label.get_text()
            words = text.split()
            if len(words) > 3:
                truncated_text = " ".join(words[:3] + ["..."])
            else:
                truncated_text = " ".join(words)
            labels.append(truncated_text)
        heatmap.set_yticklabels(labels)
        plt.tight_layout()
        
        heatmap.tick_params(left=False, bottom=False)

        # Saving graph
        figure_filename = f"{contrast}_heatmap.png"
        figure_outpath = outdir / figure_filename
        fig.savefig(figure_outpath, bbox_inches="tight")
        plt.close(fig)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-c", "--count_contrasts", type=str, required=True)
    parser.add_argument("-o", "--outdir", type=str, required=True)

    args = parser.parse_args()
    
    gsea = GeneSetEnrichmentAnalysis(Path(args.count_contrasts), Path(args.outdir))
    gsea.run()