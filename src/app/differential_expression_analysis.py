from argparse import ArgumentParser
from csv import reader
import os
import pandas as pd # type: ignore
from pathlib import Path
import pickle
from pydeseq2.dds import DeseqDataSet # type: ignore

os.environ["NUMBA_CACHE_DIR"] = "/tmp/" # Needed for scanpy to import properly
import scanpy as sc # type: ignore

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
        if samples_to_combine_path:
            gene_counts = self.combine_samples(gene_counts, samples_to_combine_path)
        gene_counts = self.rename_gene_count_ids(gene_counts, count_id_conversion)
        gene_counts = self.set_gene_id_as_index(gene_counts)
        gene_counts = self.remove_all_zero_rows(gene_counts)
        gene_counts = self.transpose_gene_counts(gene_counts)
        self.gene_counts = self.sort_by_index(gene_counts)

    @staticmethod
    def extract_gene_counts(gene_counts_path: Path) -> pd.DataFrame:
        return pd.read_csv(gene_counts_path)
    
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
    def remove_all_zero_rows(gene_counts: pd.DataFrame) -> pd.DataFrame:
        return gene_counts[gene_counts.sum(axis = 1) > 0]
    
    @staticmethod
    def transpose_gene_counts(gene_counts: pd.DataFrame) -> pd.DataFrame:
        return gene_counts.T

    @staticmethod
    def sort_by_index(gene_counts: pd.DataFrame) -> pd.DataFrame:
        return gene_counts.sort_index()

class DifferentialExpressionAnalysis:
    def __init__(self):
        sc.settings.figdir = "/src/data/pydeseq2/"
    
    def run(self, gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame) -> None:
        design_factors = ["Virus", "Time", "Lib. Prep Batch"] # Doesn't like Virus AND Treatment (redundancy probably)
        original_gene_counts, original_sample_metadata = self.remove_discrepant_sample_ids(gene_counts, sample_metadata)
        cell_lines = self.extract_cell_lines(original_sample_metadata)
        print(cell_lines)

        for cell_line in cell_lines:
            print(f"Starting on cell line: {cell_line}")
            gene_counts, sample_metadata = self.filter_for_cell_line(original_gene_counts, original_sample_metadata, cell_line)
            gene_counts, sample_metadata = self.pickle_and_rectify_batch_effect(gene_counts, sample_metadata, cell_line)

            dds = DeseqDataSet(counts=gene_counts,
                               metadata=sample_metadata,
                               design_factors=design_factors)
            dds.deseq2()
            self.perform_principal_component_analysis(dds, cell_line, design_factors)

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
    def perform_principal_component_analysis(dds: DeseqDataSet, cell_line: str, design_factors: list[str]) -> None:
        sc.tl.pca(dds)
        parameters = {"size": 200, "annotate_var_explained": True}
        for design_factor in design_factors:
            sc.pl.pca(dds, color=design_factor, save=f"_{cell_line}_{design_factor}.png", **parameters)

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