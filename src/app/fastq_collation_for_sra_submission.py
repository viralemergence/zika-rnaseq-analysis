from argparse import ArgumentParser
from csv import DictReader
from pathlib import Path
from shutil import copyfile

class BaseSpaceFastqManager:
    def __init__(self, input_dir: Path, collated_dir: Path, metadata_path: Path) -> None:
        self.input_dir = input_dir
        self.collated_dir = collated_dir
        self.metadata_path = metadata_path
        
    def run(self) -> None:
        fastq_files = self.get_basespace_fastqs(self.input_dir)
        fastq_files = self.remove_irrelevant_fastqs(fastq_files, self.metadata_path)
        
        for fastq_file in fastq_files:
            output_file = self.collated_dir / fastq_file.name
            copyfile(fastq_file, output_file)

    @staticmethod
    def get_basespace_fastqs(input_dir: Path) -> list[Path]:
        all_fastq_files = []

        subdirectories = sorted([subdir for subdir in input_dir.iterdir() if subdir.is_dir()])
        for subdir in subdirectories:
            fastq_files = [fastq for fastq in subdir.iterdir()]
            
            if len(fastq_files) > 2:
                raise Exception("More than 2 files detected in subdirectory of interest")

            for fastq in fastq_files:
                if fastq.stem.endswith("_R1_001.fastq"):
                    all_fastq_files.append(fastq)
                if fastq.stem.endswith("_R2_001.fastq"):
                    all_fastq_files.append(fastq)
        return all_fastq_files

    @classmethod
    def remove_irrelevant_fastqs(cls, fastq_files: list[Path], metadata_path: Path) -> list[Path]:
        relevant_file_names = cls.extract_relevant_file_names(metadata_path)
        return [file for file in fastq_files if file.name in relevant_file_names]

    @staticmethod
    def extract_relevant_file_names(metadata_path: Path) -> set[str]:
        file_names = set()
        with metadata_path.open() as inhandle:
            reader = DictReader(inhandle)
            for data in reader:
                if data["Cell Line"] == "HypNi":
                    continue
                r1_file_name = data["Raw_file_R1"]
                r2_file_name = r1_file_name.replace("_R1_001.fastq", "_R2_001.fastq")
                file_names.add(r1_file_name)
                file_names.add(r2_file_name)
        return file_names

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input_dir", type=str, required=True)
    parser.add_argument("-o", "--collated_dir", type=str, required=True)
    parser.add_argument("-m", "--metadata", type=str, required=True)
    args = parser.parse_args()

    bsfm = BaseSpaceFastqManager(Path(args.input_dir), Path(args.collated_dir), Path(args.metadata))
    bsfm.run()