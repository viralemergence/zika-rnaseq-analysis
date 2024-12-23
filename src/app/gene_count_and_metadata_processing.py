from argparse import ArgumentParser
from collections import defaultdict
from csv import DictWriter, reader, writer
import pandas as pd # type: ignore
from pathlib import Path
from typing import Iterator

class StarGeneCountCollater:
    def __init__(self, input_dir: Path, output_file: Path) -> None:
        self.input_files = self.get_input_files(input_dir)
        self.output_file = output_file

    @staticmethod
    def get_input_files(input_dir: Path) -> list[Path]:
        files = sorted([file for file in input_dir.iterdir() if file.is_file()])
        return files

    def run(self) -> None:
        collated_gene_counts = self.get_collated_gene_counts(self.input_files)
        gene_name_union = self.get_gene_name_union(collated_gene_counts)
        self.fill_missing_zeros(collated_gene_counts, gene_name_union)
        output_keys = ["sample_name"] + gene_name_union
        self.write_collated_gene_counts(collated_gene_counts, output_keys, self.output_file)

    @classmethod
    def get_collated_gene_counts(cls, input_files: list[Path]) -> list[defaultdict[int]]:
        collated_gene_counts = list()
        for file in input_files:
            sample_info = defaultdict(int)
            sample_info["sample_name"] = cls.extract_sample_name(file)
            sample_info.update(cls.extract_gene_counts(file))
            collated_gene_counts.append(sample_info)
        return collated_gene_counts
        
    @staticmethod
    def extract_sample_name(file: Path) -> str:
        suffix = "_L001_ReadsPerGene.out.tab"
        return file.name.replace(suffix, "")
    
    @classmethod
    def extract_gene_counts(cls, file: Path) -> dict[int]:
        gene_counts = dict()
        with file.open() as inhandle:
            reader_iterator = reader(inhandle, delimiter="\t")
            cls.skip_gene_count_header(reader_iterator)
            for line in reader_iterator:
                gene = line[0]
                strand_1_count = int(line[2])
                strand_2_count = int(line[3])
                total_count = strand_1_count + strand_2_count

                gene_counts[gene] = total_count
        return gene_counts

    @staticmethod
    def skip_gene_count_header(reader_iterator: Iterator) -> None:
        header_set = {"N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"}
        while len(header_set) > 0:
            header = next(reader_iterator)[0]
            header_set.remove(header)
            
    @staticmethod
    def get_gene_name_union(collated_gene_counts: list[defaultdict[int]]) -> list[str]:
        gene_names = set()
        for sample in collated_gene_counts:
            for gene, _ in sample.items():
                gene_names.add(gene)
        gene_names.remove("sample_name")
        return sorted(list(gene_names))

    @staticmethod
    def fill_missing_zeros(collated_gene_counts: list[defaultdict[int]], gene_name_union: list[str]) -> None:
        for sample in collated_gene_counts:
            for gene in gene_name_union:
                sample[gene] = sample[gene]
    
    @staticmethod
    def write_collated_gene_counts(collated_gene_counts: list[defaultdict[int]], output_keys: list[str], output_file: Path) -> None:
        with output_file.open("w") as outhandle:
            writer = DictWriter(outhandle, output_keys)
            writer.writeheader()
            writer.writerows(collated_gene_counts)
            
    def transpose_gene_counts_table(self) -> None:
        with self.output_file.open() as inhandle:
            data = [line for line in reader(inhandle)]
            data = self.transpose_list_of_lists(data)

        with self.output_file.open("w") as outhandle:
            data_writer = writer(outhandle)
            data_writer.writerows(data)
            
    @staticmethod
    def transpose_list_of_lists(data: list[list]) -> list[list]:
        return [list(i) for i in zip(*data)]

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
        gene_counts = self.remove_blacklist_samples(gene_counts)
        if samples_to_combine_path:
            gene_counts = self.combine_samples(gene_counts, samples_to_combine_path)
        gene_counts = self.rename_gene_count_ids(gene_counts, count_id_conversion)
        gene_counts = self.set_gene_id_as_index(gene_counts)
        gene_counts = self.remove_low_count_rows(gene_counts)
        gene_counts = self.transpose_gene_counts(gene_counts)
        self.gene_counts = self.sort_by_index(gene_counts)

    @staticmethod
    def extract_gene_counts(gene_counts_path: Path) -> pd.DataFrame:
        return pd.read_csv(gene_counts_path)
    
    @staticmethod
    def remove_blacklist_samples(gene_counts: pd.DataFrame) -> pd.DataFrame:
        blacklist_samples = ["HypNi_ZIKV_PRVABC59_24_a_S35", "Julianna-87_S31"]
        return gene_counts.drop(blacklist_samples, axis=1)

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
    def remove_low_count_rows(gene_counts: pd.DataFrame, threshold: int = 10) -> pd.DataFrame:
        return gene_counts[gene_counts.sum(axis = 1) > threshold]
    
    @staticmethod
    def transpose_gene_counts(gene_counts: pd.DataFrame) -> pd.DataFrame:
        return gene_counts.T

    @staticmethod
    def sort_by_index(gene_counts: pd.DataFrame) -> pd.DataFrame:
        return gene_counts.sort_index()

class FinalFormatter:
    def __init__(self) -> None:
        pass

    def run(self, gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame, out_dir: Path) -> None:
        gene_counts, sample_metadata = self.remove_discrepant_sample_ids(gene_counts, sample_metadata)
        self.write_table(gene_counts, out_dir, "processed_gene_counts.csv")
        self.write_table(sample_metadata, out_dir, "processed_sample_metadata.csv")

    @staticmethod
    def remove_discrepant_sample_ids(gene_counts: pd.DataFrame, sample_metadata: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
        gene_count_sample_ids = set(gene_counts.index)
        metadata_sample_ids = set(sample_metadata.index)
        shared_ids = gene_count_sample_ids & metadata_sample_ids

        gene_counts = gene_counts[gene_counts.index.isin(shared_ids)]
        sample_metadata = sample_metadata[sample_metadata.index.isin(shared_ids)]
        return gene_counts, sample_metadata
    
    @staticmethod
    def write_table(table: pd.DataFrame, out_dir: Path, outname: str) -> None:
        outpath = out_dir / outname
        table.to_csv(outpath)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input_dir", type=str, required=True)
    parser.add_argument("-c", "--collated_gene_counts", type=str, required=True)
    parser.add_argument("-m", "--metadata", type=str, required=True)
    parser.add_argument("-k", "--combine", type=str, required=False)
    args = parser.parse_args()

    sgcc = StarGeneCountCollater(Path(args.input_dir), Path(args.collated_gene_counts))
    sgcc.run()
    sgcc.transpose_gene_counts_table()

    smdm = SampleMetaDataManager(Path(args.metadata))
    smdm.run()

    gcm = GeneCountManager(Path(args.collated_gene_counts))
    gcm_parameters = {"count_id_conversion": smdm.count_id_conversion}
    if args.combine:
        gcm_parameters.update({"samples_to_combine_path": Path(args.combine)})
    gcm.run(**gcm_parameters)
    
    ff = FinalFormatter()
    ff.run(gcm.gene_counts, smdm.sample_metadata, gcm.gene_counts_path.parent)