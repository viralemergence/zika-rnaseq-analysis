from argparse import ArgumentParser
from csv import reader
import matplotlib.pyplot as plt # type: ignore
from matplotlib_venn import venn2 # type: ignore
from pathlib import Path

class LRTVennDiagram:
    def __init__(self, cell_lines: list[str], viruses: list[str], lrt_dir: Path, outdir: Path) -> None:
        self.cell_lines = cell_lines
        self.viruses = viruses
        self.lrt_dir = lrt_dir
        self.outdir = outdir

    def run(self) -> None:
        gene_sets = {}
        for cell_line in self.cell_lines:
            print(f"\nStarting on cell line: {cell_line}")
            for virus in self.viruses:
                print(f"Starting on virus: {virus}")
                
                lrt_results = self.lrt_dir / f"{cell_line}_{virus}_LRT_results.csv"
                gene_set = self.extract_significant_lrt_genes(lrt_results)
                gene_sets[f"{cell_line}_{virus}"] = gene_set

        set1 = gene_sets["R06E_MR"]
        set2 = gene_sets["R06E_PRV"]
        
        fig, ax = plt.subplots()
        
        venn_diagram = venn2([set1, set2], set_labels=["R06E MR", "R06E PRV"])
        
        circle_coordinates = (0.2, -0.75)
        circle_test = plt.Circle(circle_coordinates, 0.1, facecolor="lightblue", edgecolor="none", clip_on=False)
        ax.add_patch(circle_test)
        ax.text(circle_coordinates[0], circle_coordinates[1], "0", ha="center", va="center")
        ax.text(circle_coordinates[0], circle_coordinates[1]-.15, "Aji MR & PRV", ha="center", va="center")

        figure_filename = f"lrt_results_venn_diagram.png"
        figure_outpath = self.outdir / figure_filename
        fig.savefig(figure_outpath, bbox_inches="tight")
        plt.close(fig)
    
    @staticmethod
    def extract_significant_lrt_genes(lrt_results: Path) -> set[str]:
        gene_set = set()
        with lrt_results.open() as inhandle:
            reader_iterator = reader(inhandle, delimiter=",")
            header = next(reader_iterator)
            for line in reader_iterator:
                try:
                    padj = float(line[-1])
                except ValueError:
                    continue
                if padj > 0.05:
                    continue
                gene_set.add(line[0])
        return gene_set

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-c", "--cell_lines", nargs="+", required=True)
    parser.add_argument("-v", "--viruses", nargs="+", required=True)
    parser.add_argument("-l", "--lrt_dir", type=str, required=True)
    parser.add_argument("-o", "--outdir", type=str, required=True)

    args = parser.parse_args()
    
    lvd = LRTVennDiagram(args.cell_lines, args.viruses, Path(args.lrt_dir), Path(args.outdir))
    lvd.run()