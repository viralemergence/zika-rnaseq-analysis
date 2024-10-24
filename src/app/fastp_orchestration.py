from argparse import ArgumentParser
from pathlib import Path

class FastqInputManager:
    def __init__(self, input_dir: Path, output_dir: Path, subdir_index: int) -> None:
        self.input_fastqs = self.get_input_fastqs(input_dir, subdir_index)
        self.output_fastqs = self.set_output_fastqs(self.input_fastqs, output_dir)
        
    @staticmethod
    def get_input_fastqs(input_dir: Path, subdir_index: int) -> dict[Path]:
        subdirectories = sorted([subdir for subdir in input_dir.iterdir() if subdir.is_dir()])
        subdir = subdirectories[subdir_index]
        fastq_files = [fastq for fastq in subdir.iterdir()]
        
        if len(fastq_files) > 2:
            raise Exception("More than 2 files detected in subdirectory of interest")
        
        input_fastqs = dict()
        for fastq in fastq_files:
            if fastq.stem.endswith("_R1_001.fastq"):
                input_fastqs["R1"] = fastq
            if fastq.stem.endswith("_R2_001.fastq"):
                input_fastqs["R2"] = fastq
        return input_fastqs
        
    @staticmethod
    def set_output_fastqs(input_fastqs: dict[Path], output_dir: Path) -> dict[Path]:
        output_fastqs = dict()
        for key, path in input_fastqs.items():
            outpath_stem = path.stem.replace(".fastq", "")
            outpath = output_dir / f"{outpath_stem}_fastp.fastq.gz"
            output_fastqs[key] = outpath
        return output_fastqs

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input_dir", type=str, required=True)
    parser.add_argument("-o", "--output_dir", type=str, required=True)
    parser.add_argument("-j", "--subdir_index", type=int, required=True)
    args = parser.parse_args()

    fim = FastqInputManager(Path(args.input_dir), Path(args.output_dir), args.subdir_index)
    print()
    print(fim.input_fastqs)
    print(fim.output_fastqs)