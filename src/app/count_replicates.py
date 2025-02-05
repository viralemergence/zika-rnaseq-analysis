from argparse import ArgumentParser
from collections import defaultdict
from csv import DictReader
from pathlib import Path

class ReplicateCounter:
    def __init__(self, metadata_path: Path, outdir: Path) -> None:
        self.metadata_path = metadata_path
        self.outdir = outdir

    def run(self) -> None:
        replicate_counts = self.calculate_replicate_counts(self.metadata_path)

        outpath = self.outdir / f"replicate_counts.csv"
        with outpath.open("w") as outhandle:
            for key, value in replicate_counts.items():
                outline = ",".join([str(s) for s in [key, value]])
                outhandle.write(outline + "\n")

    @staticmethod
    def calculate_replicate_counts(metadata_path: Path) -> defaultdict[int]:
        replicate_counts = defaultdict(int)
        with metadata_path.open() as inhandle:
            metadata = DictReader(inhandle, delimiter=",")
            for sample in metadata:
                cell_line = sample["Cell Line"]
                if cell_line == "HypNi":
                    continue
                virus = sample["Virus"]
                time = sample["Time"]
                if time == "0.0":
                    continue
                sample_type = "-".join([cell_line, virus, time])
                replicate_counts[sample_type] += 1
        return replicate_counts

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-m", "--metadata", type=str, required=True)
    parser.add_argument("-o", "--outdir", type=str, required=True)

    args = parser.parse_args()
    
    rc = ReplicateCounter(Path(args.metadata), Path(args.outdir))
    rc.run()