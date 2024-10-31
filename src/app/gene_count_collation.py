from argparse import ArgumentParser
from collections import defaultdict
from csv import DictWriter, reader, writer
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

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input_dir", type=str, required=True)
    parser.add_argument("-o", "--output_file", type=str, required=True)
    args = parser.parse_args()

    sgcc = StarGeneCountCollater(Path(args.input_dir), Path(args.output_file))
    sgcc.run()
    sgcc.transpose_gene_counts_table()
