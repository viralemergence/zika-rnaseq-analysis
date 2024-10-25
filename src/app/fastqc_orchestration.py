from argparse import ArgumentParser
from pathlib import Path
import subprocess

class FastqPathManager:
    def __init__(self, input_dir: Path, file_index: int) -> None:
        self.input_fastq = self.set_input_fastq(input_dir, file_index)
        
    @staticmethod
    def set_input_fastq(input_dir: Path, file_index: int) -> Path:
        fastq_files = sorted([file for file in input_dir.iterdir() if file.is_file()])
        return fastq_files[file_index]

class FastqcManager:
    def __init__(self, input_fastq: Path, output_dir: Path) -> None:
        self.input_fastq = input_fastq
        self.output_dir = output_dir

    def run_fastqc(self) -> None:
        fastqc_command = ["fastqc", "-o", self.output_dir, self.input_fastq]

        print(self.input_fastq)
        p = subprocess.Popen(fastqc_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        while p.poll() is None and (line := p.stderr.readline()) != "":
            print(line.strip())
        p.wait()
        print(f"Exit code: {p.poll()}")

        if p.poll() != 0:
            raise Exception("Fastqc did not complete successfully")

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input_dir", type=str, required=True)
    parser.add_argument("-o", "--output_dir", type=str, required=True)
    parser.add_argument("-j", "--file_index", type=int, required=True)
    args = parser.parse_args()

    fpm = FastqPathManager(Path(args.input_dir), args.file_index)
    fpm = FastqcManager(fpm.input_fastq, Path(args.output_dir))
    fpm.run_fastqc()