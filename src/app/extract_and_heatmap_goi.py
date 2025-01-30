from argparse import ArgumentParser
from collections import defaultdict
from csv import reader
from numpy import log10, where # type: ignore
import pandas as pd # type: ignore
from pathlib import Path

class ContrastGOIManager:
    def __init__(self, count_contrasts_path: Path, goi_path: Path, cell_line: str, outdir: Path, virus_contrast: str) -> None:
        self.count_contrasts_path = count_contrasts_path
        self.genes_of_interest = self.extract_genes_of_interest(goi_path)
        self.cell_line = cell_line
        self.outdir = outdir
        self.virus_contrast = virus_contrast

        self.check_count_contrasts_filename(self.count_contrasts_path, self.cell_line)

    @staticmethod
    def extract_genes_of_interest(goi_path: Path) -> list[str]:
        with goi_path.open() as inhandle:
            reader_iterator = reader(inhandle, delimiter=",")
            return list(map(lambda line: line[0], reader_iterator))

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
        goi_fold_changes = self.extract_goi_log2_fold_change(count_contrasts, self.genes_of_interest, self.virus_contrast)
        goi_fold_changes = self.load_significant_genes_to_pd_dataframe(goi_fold_changes)
        print(goi_fold_changes)
        quit()

        significant_genes_outpath = self.outdir / f"{self.cell_line}_deg_genes_by_contrast.csv"
        significant_genes.to_csv(significant_genes_outpath, index=True)

    @staticmethod
    def extract_count_contrasts(count_contrasts: Path) -> pd.DataFrame:
        df = pd.read_csv(count_contrasts, index_col=0)
        columns = [eval(column) for column in df.columns]
        df.columns = pd.MultiIndex.from_tuples(columns)
        return df

    @staticmethod
    def extract_goi_log2_fold_change(count_contrasts: pd.DataFrame, genes_of_interest: list[str], virus_contrast: str) -> dict[pd.DataFrame]:
        log2_fold_changes = dict()
        for contrast, deseq_results in count_contrasts.T.groupby(level=0):
            if contrast != virus_contrast:
                continue
            deseq_results = deseq_results.T.droplevel(0, axis=1)
            deseq_results.columns = deseq_results.columns.set_levels(deseq_results.columns.levels[0].astype(float), level=0)

            for time, deseq_results_by_time in deseq_results.T.groupby(level=0):
                deseq_results_by_time = deseq_results_by_time.T.droplevel(0, axis=1)

                log2_fold_change = deseq_results_by_time["log2FoldChange"].loc[genes_of_interest]
                log2_fold_changes[time] = log2_fold_change
        return log2_fold_changes
    
    @staticmethod
    def load_significant_genes_to_pd_dataframe(goi_fold_changes: dict[pd.DataFrame]) -> pd.DataFrame:
        dataframe_list = []
        columns = []
        for time_point, data in goi_fold_changes.items():
            columns.append(time_point)
            dataframe_list.append(data)
        df = pd.concat(dataframe_list, axis=1).sort_index()
        df.columns = columns
        return df

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-k", "--count_contrasts", type=str, required=True)
    parser.add_argument("-g", "--goi", type=str, required=True)
    parser.add_argument("-l", "--cell_line", type=str, required=True)
    parser.add_argument("-o", "--outdir", type=str, required=True)
    parser.add_argument("-v", "--virus_contrast", type=str, required=True)

    args = parser.parse_args()
    
    cdm = ContrastGOIManager(Path(args.count_contrasts), Path(args.goi), args.cell_line, Path(args.outdir), args.virus_contrast)
    cdm.run()