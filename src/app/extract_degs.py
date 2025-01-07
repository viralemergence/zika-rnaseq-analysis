from argparse import ArgumentParser
from collections import defaultdict
from numpy import log10, where # type: ignore
import pandas as pd # type: ignore
from pathlib import Path

class ContrastDEGManager:
    def __init__(self, count_contrasts_path: Path, cell_line: str, outdir: Path) -> None:
        self.count_contrasts_path = count_contrasts_path
        self.cell_line = cell_line
        self.outdir = outdir

        self.check_count_contrasts_filename(self.count_contrasts_path, self.cell_line)

    @staticmethod
    def check_count_contrasts_filename(count_contrasts_path: Path, cell_line: str) -> None:
        if cell_line in count_contrasts_path.stem:
            print(f"\nCell line '{cell_line}' detected in file name.\
                  \nProceeding with extraction")
            return
        else:
            print(f"\nWARNING: Cell line '{cell_line}' not detected in file name.\
                  \nExtraction continuing, but double check correct input is being used")

    def run(self) -> None:
        count_contrasts = self.extract_count_contrasts(self.count_contrasts_path)
        significant_genes = self.extract_significant_genes(count_contrasts)
        significant_genes = self.load_significant_genes_to_pd_dataframe(significant_genes)

        significant_genes_outpath = self.outdir / f"{self.cell_line}_deg_genes_by_contrast.csv"
        significant_genes.to_csv(significant_genes_outpath, index=True)

    @staticmethod
    def extract_count_contrasts(count_contrasts: Path) -> pd.DataFrame:
        df = pd.read_csv(count_contrasts, index_col=0)
        columns = [eval(column) for column in df.columns]
        df.columns = pd.MultiIndex.from_tuples(columns)
        return df

    @staticmethod
    def extract_significant_genes(count_contrasts: pd.DataFrame) -> defaultdict[dict[pd.DataFrame]]:
        contrasts_significant_genes = defaultdict(dict)
        for contrast, deseq_results in count_contrasts.T.groupby(level=0):
            deseq_results = deseq_results.T.droplevel(0, axis=1)
            deseq_results.columns = deseq_results.columns.set_levels(deseq_results.columns.levels[0].astype(float), level=0)

            for time, deseq_results_by_time in deseq_results.T.groupby(level=0):
                deseq_results_by_time = deseq_results_by_time.T.droplevel(0, axis=1)
                deseq_results_by_time["Significant"] = where((abs(deseq_results_by_time["log2FoldChange"]) > 1.5) & (deseq_results_by_time["padj"] < 0.05), "yes", "no")
                significant_genes = deseq_results_by_time.loc[deseq_results_by_time["Significant"] == "yes"].copy()
                significant_genes["Rank"] = -log10(significant_genes["padj"]*abs(significant_genes["log2FoldChange"]))
                significant_genes = significant_genes.drop(columns=["Significant"])
                
                contrasts_significant_genes[contrast][time] = significant_genes
        return contrasts_significant_genes
    
    @staticmethod
    def load_significant_genes_to_pd_dataframe(significant_genes: defaultdict[dict[pd.DataFrame]]) -> pd.DataFrame:
        dataframe_list = []
        for contrast, time in significant_genes.items():
            for time_point, data in time.items():
                new_column_names = [f"({contrast}, {time_point}, {column})" for column in data.columns]
                data.columns = new_column_names
                dataframe_list.append(data)
        return pd.concat(dataframe_list, axis=1).sort_index()

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-k", "--count_contrasts", type=str, required=True)
    parser.add_argument("-l", "--cell_line", type=str, required=True)
    parser.add_argument("-o", "--outdir", type=str, required=True)

    args = parser.parse_args()
    
    cdm = ContrastDEGManager(Path(args.count_contrasts), args.cell_line, Path(args.outdir))
    cdm.run()