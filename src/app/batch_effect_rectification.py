from argparse import ArgumentParser
from inmoose.pycombat import pycombat_seq # type: ignore
import pandas as pd # type: ignore
from pathlib import Path

class BatchEffectRectifier:
    def run(self, gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame, out_dir: Path) -> None:
        cell_lines = self.extract_cell_lines(sample_metadata)
        print(f"Metadata listed cell lines: {cell_lines}")

        gene_counts_by_cell_line_rectified = []
        sample_metadata_by_cell_line_rectified = []

        for cell_line in cell_lines:
            print(f"\nStarting on cell line: {cell_line}\n----------")
            gene_counts_by_cell_line, sample_metadata_by_cell_line = self.filter_for_cell_line(gene_counts, sample_metadata, cell_line)
            gene_counts_by_cell_line, sample_metadata_by_cell_line = self.rectify_batch_effect(gene_counts_by_cell_line, sample_metadata_by_cell_line)
            
            gene_counts_by_cell_line_rectified.append(gene_counts_by_cell_line)
            sample_metadata_by_cell_line_rectified.append(sample_metadata_by_cell_line)
            
        gene_counts_rectified = pd.concat(gene_counts_by_cell_line_rectified, axis=0)
        sample_metadata_rectified = pd.concat(sample_metadata_by_cell_line_rectified, axis=0)
        
        self.write_table(gene_counts_rectified, out_dir, "gene_counts_final.csv")
        self.write_table(sample_metadata_rectified, out_dir, "sample_metadata_final.csv")

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
    def rectify_batch_effect(cls, gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
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
    def write_table(table: pd.DataFrame, out_dir: Path, outname: str) -> None:
        outpath = out_dir / outname
        table.to_csv(outpath)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-c", "--counts", type=str, required=True)
    parser.add_argument("-m", "--metadata", type=str, required=True)
    parser.add_argument("-o", "--out_dir", type=str, required=True)
    args = parser.parse_args()

    gene_counts = pd.read_csv(args.counts, index_col=0)
    sample_metadata = pd.read_csv(args.metadata, index_col=0)
    out_dir = Path(args.out_dir)
    
    r = BatchEffectRectifier()
    r.run(gene_counts, sample_metadata, out_dir)