from argparse import ArgumentParser
from csv import reader
import matplotlib.pyplot as plt # type: ignore
from matplotlib_venn import venn3
from pathlib import Path
import seaborn as sns # type: ignore

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

        print(len(gene_sets))

        return            
        for f in "Uu":
            fig, ax = plt.subplots()
        
            heatmap = sns.heatmap(correlation_matrix, ax=ax, vmin=0.9)
            cbar = heatmap.collections[0].colorbar
            cbar.set_label("Spearman Coefficient", labelpad=10)
            columns = [4, 8]
            for col in columns:
                heatmap.axvline(col, color="white", lw=3)
                heatmap.axhline(col, color="white", lw=3)
            heatmap.set_xlabel("")
            heatmap.set_ylabel("")
            
            heatmap.tick_params(left=False, bottom=False)

            figure_filename = f"{cell_line}_{method}_heatmap.png"
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