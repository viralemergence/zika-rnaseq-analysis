from argparse import ArgumentParser
from csv import reader
from pathlib import Path
from shutil import move
import subprocess

class FastqInputManager:
    def __init__(self, input_dir: Path, file_names: Path, file_index: int) -> None:
        self.input_fastqs = self.get_input_fastqs(input_dir, file_names, file_index)

    @classmethod
    def get_input_fastqs(cls, input_dir: Path, file_names: Path, file_index: int) -> dict[Path]:
        file_names = cls.get_file_names(file_names)
        r1_fastq_stem = file_names[file_index]
        
        input_fastqs = dict()
        input_fastqs["R1"] = input_dir / r1_fastq_stem
        input_fastqs["R2"] = input_dir / r1_fastq_stem.replace("_R1_", "_R2_")
        return input_fastqs
    
    @staticmethod
    def get_file_names(file_names: Path) -> list[str]:
        files = []
        with file_names.open() as inhandle:
            reader_iterator = reader(inhandle)
            for line in reader_iterator:
                files.append(line[0])
        return files

class StarManager:
    def __init__(self, input_fastqs: dict[Path], genome_index: Path, output_dir: Path) -> None:
        self.input_fastqs = input_fastqs
        self.genome_index = genome_index
        self.sub_output_dir = self.set_sub_output_dir(output_dir, input_fastqs["R1"])

    @staticmethod
    def set_sub_output_dir(output_dir: Path, r1_path: Path) -> Path:
        sub_output_dir = output_dir / (r1_path.stem.replace("_R1_001_fastp.fastq", "")+"/")
        sub_output_dir.mkdir(exist_ok=True)
        return sub_output_dir

    def run_star_gene_count(self) -> None:
        input_r1 = self.input_fastqs["R1"]
        input_r2 = self.input_fastqs["R2"]

        star_command = ["STAR",
                         "--genomeDir", self.genome_index, "--runThreadN", "10",
                         "--readFilesCommand", "zcat",
                         "--readFilesIn", input_r1, input_r2,
                         "--outFileNamePrefix", f"{self.sub_output_dir}/",
                         "--quantMode", "GeneCounts"
                         ]

        p = subprocess.Popen(star_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        while p.poll() is None and (line := p.stdout.readline()) != "":
            print(line.strip())
        p.wait()
        print(f"Exit code: {p.poll()}")

        if p.poll() != 0:
            print(p.stderr.readlines())
            raise Exception("STAR did not complete successfully")
        
    def move_gene_counts(self) -> None:
        gene_count_path = self.sub_output_dir / "ReadsPerGene.out.tab"
        new_file_path = "/src/data/star/gene_counts/" + self.sub_output_dir.name + "_ReadsPerGene.out.tab"
        move(gene_count_path, new_file_path)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-i", "--input_dir", type=str, required=True)
    parser.add_argument("-f", "--file_names", type=str, required=True)
    parser.add_argument("-g", "--genome_index", type=str, required=True)
    parser.add_argument("-o", "--output_dir", type=str, required=True)
    parser.add_argument("-j", "--file_index", type=int, required=True)
    args = parser.parse_args()

    fim = FastqInputManager(Path(args.input_dir), Path(args.file_names), args.file_index)
    sm = StarManager(fim.input_fastqs, Path(args.genome_index), Path(args.output_dir))
    sm.run_star_gene_count()
    sm.move_gene_counts()