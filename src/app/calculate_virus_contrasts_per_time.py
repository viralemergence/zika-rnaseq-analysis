from argparse import ArgumentParser
from contextlib import contextmanager, redirect_stdout, redirect_stderr
from os import devnull
import pandas as pd # type: ignore
from pathlib import Path
from pydeseq2.ds import DeseqStats # type: ignore
from pydeseq2.dds import DeseqDataSet # type: ignore
from warnings import catch_warnings, simplefilter

@contextmanager
def suppress_stdout_stderr():
    with open(devnull, "w") as null_handle:
        with redirect_stdout(null_handle) as out, redirect_stderr(null_handle) as err:
            yield(out, err)

class CalculateContrasts:
    def __init__(self, cell_lines: list[str], virus_contrasts: list[str], outdir: Path) -> None:
        self.cell_lines = cell_lines
        self.virus_contrasts = [contrast.split("vs") for contrast in virus_contrasts]
        self.outdir = outdir

    def run(self, gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame) -> None:
        time_points_of_interest = [6.0, 12.0, 24.0, 48.0]
        for cell_line in self.cell_lines:
            print(f"\nStarting on cell line: {cell_line}")
            gene_counts_by_cell_line, sample_metadata_by_cell_line = self.filter_for_cell_line(gene_counts, sample_metadata, cell_line)
            
            contrast_results = []
            contrast_labels = []
            
            for treatment, control in self.virus_contrasts:
                print(f"\nStarting on virus contrast: {treatment} vs {control}")
                gene_counts_by_virus, sample_metadata_by_virus = self.filter_for_virus_contrasts(gene_counts_by_cell_line, sample_metadata_by_cell_line, treatment, control)
                
                for time_point in time_points_of_interest:
                    print(f"Starting on time: {time_point}")
                    gene_counts_by_time, sample_metadata_by_time = self.filter_for_time(gene_counts_by_virus, sample_metadata_by_virus, time_point)

                    dds = self.pydeseq2_normalization(gene_counts_by_time, sample_metadata_by_time, ["Virus"])
                    deseq_results = self.deseq_stats_calculation(dds, treatment, control)
                    self.extract_and_append_foldchange_padj(deseq_results, contrast_results, contrast_labels, treatment, control, time_point)

            collated_df = pd.concat(contrast_results, axis=1)
            collated_df.columns = contrast_labels
            print("\nSaving collated results")
            print(collated_df)
            
            out_filename = f"{cell_line}_virus_contrasts_per_time.csv"
            outpath = self.outdir / out_filename
            collated_df.to_csv(outpath)
    
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
    def filter_for_time(gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame, time: float) -> tuple[pd.DataFrame, pd.DataFrame]:
        sample_metadata = sample_metadata[sample_metadata["Time"] == time]
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
    def deseq_stats_calculation(dds: DeseqDataSet, treatment: str, control: str) -> pd.DataFrame:
        results = DeseqStats(dds, contrast=["Virus", treatment, control.replace("_", "-")])
        with suppress_stdout_stderr():
            results.summary()
        return results.results_df

    @staticmethod
    def extract_and_append_foldchange_padj(deseq_results: pd.DataFrame, contrast_results: list[pd.DataFrame], contrast_labels: list[str],
                                           treatment: str, control: str, time_point: float) -> None:
        columns_of_interest = ["log2FoldChange", "padj"]
        for column in columns_of_interest:
            contrast_results.append(deseq_results[column])
            label = (f"{treatment}_vs_{control}", f"{time_point}", column)
            contrast_labels.append(label)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-c", "--counts", type=str, required=True)
    parser.add_argument("-m", "--metadata", type=str, required=True)
    parser.add_argument("-l", "--cell_lines", nargs="+", required=True)
    parser.add_argument("-v", "--virus_contrasts", nargs="+", required=True)
    parser.add_argument("-o", "--outdir", type=str, required=True)

    args = parser.parse_args()

    gene_counts = pd.read_csv(args.counts, index_col=0)
    sample_metadata = pd.read_csv(args.metadata, index_col=0)
    
    cc = CalculateContrasts(args.cell_lines, args.virus_contrasts, Path(args.outdir))
    cc.run(gene_counts, sample_metadata)