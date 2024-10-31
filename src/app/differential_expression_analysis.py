from argparse import ArgumentParser
import pandas as pd
from pathlib import Path

class GeneCountManager:
    def __init__(self, gene_counts_path: Path) -> None:
        self.gene_counts_path = gene_counts_path

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-c", "--counts", type=str, required=True)
    args = parser.parse_args()

    gcm = GeneCountManager(Path(args.counts))
