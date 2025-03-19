from argparse import ArgumentParser
from csv import reader
from matplotlib import cm # type: ignore
from matplotlib.colors import LinearSegmentedColormap # type: ignore
import matplotlib.pyplot as plt # type: ignore
import networkx as nx # type: ignore
import numpy as np # type: ignore
import pandas as pd # type: ignore
from pathlib import Path
import requests

class StringNetworkManager:
    def __init__(self, goi: Path, outdir: Path) -> None:
        self.genes_of_interest = self.extract_genes_of_interest(goi)
        self.outdir = outdir

    @staticmethod
    def extract_genes_of_interest(goi_path: Path) -> list[str]:
        with goi_path.open() as inhandle:
            reader_iterator = reader(inhandle, delimiter=",")
            return list(map(lambda line: line[0], reader_iterator))

    def run(self) -> None:
        goi_submission = "%0d".join(self.genes_of_interest)
        taxon_id = 9606 # Note: Using human for now
        url = f"https://string-db.org/api/tsv/network?identifiers={goi_submission}&species={taxon_id}"
        r = requests.get(url)
        lines = r.text.split("\n")
        data = [line.split("\t") for line in lines]
        df = pd.DataFrame(data[1:-1], columns=data[0])
        interactions = df[['preferredName_A', 'preferredName_B', 'score']]

        G=nx.Graph(name='Protein Interaction Graph')
        interactions = np.array(interactions)
        for i in range(len(interactions)):
            interaction = interactions[i]
            a = interaction[0]
            b = interaction[1]
            w = float(interaction[2])
            G.add_weighted_edges_from([(a,b,w)])

        mean_connection_count = (G.number_of_edges()*2) / G.number_of_nodes()
        print(mean_connection_count)

        cmap_colors = ["red", "blue"][::-1]
        cmap = LinearSegmentedColormap.from_list("red_to_purple", cmap_colors)

        def rescale(l,newmin,newmax):
            arr = list(l)
            return [(x-min(arr))/(max(arr)-min(arr))*(newmax-newmin)+newmin for x in arr]

        nc = rescale([G.degree(v) for v in G],0.0,0.9) 
        nc = [cmap(i) for i in nc]

        pos = nx.spring_layout(G, k=2)

        plt.rcParams["svg.fonttype"] = "none"
        fig, ax = plt.subplots(figsize=(6,5), dpi=1200)
        nx.draw_networkx(G, pos,
                         node_color=nc, node_shape="o", node_size=1300,
                         font_color="white", font_size=8, font_weight="bold",
                         edge_color="grey", width=1
                         )
        plt.axis('off')

        sm = cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 16))
        sm._A=[]
        cbar = plt.colorbar(sm, ax=ax, location="right", shrink=0.5, pad=0)
        cbar.ax.set_title("PPI\nCount")#, fontweight="bold")

        plt.tight_layout()
        figure_outpath = self.outdir / "ppi_network.png"
        plt.savefig(figure_outpath, transparent=True)
        figure_outpath = self.outdir / "ppi_network.svg"
        plt.savefig(figure_outpath, transparent=True)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-g", "--goi", type=str, required=True)
    parser.add_argument("-o", "--outdir", type=str, required=True)

    args = parser.parse_args()
    
    snm = StringNetworkManager(Path(args.goi), Path(args.outdir))
    snm.run()