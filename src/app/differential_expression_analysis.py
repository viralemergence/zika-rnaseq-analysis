from argparse import ArgumentParser
import pandas as pd
from pathlib import Path
from pydeseq2.dds import DeseqDataSet

class SampleMetaDataManager:
    def __init__(self, metadata_path: Path) -> None:
        self.metadata_path = metadata_path
        
    def run(self) -> None:
        sample_metadata = self.extract_sample_metadata(self.metadata_path)
        sample_metadata = self.remove_rows_without_id(sample_metadata)
        sample_metadata = self.set_gene_count_id(sample_metadata)
        sample_metadata = self.set_sample_id_as_index(sample_metadata)
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
    def set_count_id_to_sample_id_dict(sample_metadata: pd.DataFrame) -> dict[str]:
        gene_count_ids = sample_metadata["gene_count_id"].tolist()
        sample_ids = sample_metadata.index.tolist()
        return {count_id: sample_id for count_id, sample_id in zip(gene_count_ids, sample_ids)}

class GeneCountManager:
    def __init__(self, gene_counts_path: Path) -> None:
        self.gene_counts_path = gene_counts_path

    def run(self, count_id_conversion: dict[str]) -> None:
        gene_counts = self.extract_gene_counts(self.gene_counts_path)
        gene_counts = self.rename_gene_count_ids(gene_counts, count_id_conversion)
        gene_counts = self.set_gene_id_as_index(gene_counts)
        gene_counts = self.remove_all_zero_rows(gene_counts)
        self.gene_counts = self.transpose_gene_counts(gene_counts)

    @staticmethod
    def extract_gene_counts(gene_counts_path: Path) -> pd.DataFrame:
        return pd.read_csv(gene_counts_path)

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

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-c", "--counts", type=str, required=True)
    parser.add_argument("-m", "--metadata", type=str, required=True)
    args = parser.parse_args()

    smdm = SampleMetaDataManager(Path(args.metadata))
    smdm.run()

    gcm = GeneCountManager(Path(args.counts))
    gcm.run(smdm.count_id_conversion)