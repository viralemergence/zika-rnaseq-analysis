from collections import defaultdict
from csv import reader
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.base import download_go_basic_obo, dnld_file
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.obo_parser import GODag
from pathlib import Path

class GoatoolsManager:
    def __init__(self) -> None:
        self.background_genes_path = Path("/src/data/pydeseq2/goatools/ncbi_gene_results.txt")
        self.obo_path = Path("/src/data/pydeseq2/goatools/go-basic.obo")
        self.gene2go_ftp_path = "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz"
        self.gene2go_path = Path("/src/data/pydeseq2/goatools/gene2go.txt")
        self.study_genes_path = Path("/src/data/pydeseq2/degpatterns/gene_clusters.csv")

    def setup(self) -> None:
        download_go_basic_obo(str(self.obo_path))
        dnld_file(self.gene2go_ftp_path, str(self.gene2go_path))
        self.get_background_genes(self.background_genes_path)
        
    def run(self) -> None:
        background_genes = self.extract_background_genes(self.background_genes_path)
        symbol_id_mapper = self.set_symbol_id_mapper(self.background_genes_path)
        study_genes = self.extract_study_genes(self.study_genes_path)
        study_genes = self.convert_study_gene_symbols(study_genes, symbol_id_mapper)

        godag = GODag(self.obo_path)

        objanno = Gene2GoReader(str(self.gene2go_path), taxids=[9407])
        ns2assoc = objanno.get_ns2assc()
        id2gos = defaultdict(set)
        for _, associations in ns2assoc.items():
            for gene_id, go_terms in associations.items():
                id2gos[gene_id].update(go_terms)

        goeaobj = GOEnrichmentStudy(
            background_genes,
            id2gos,
            godag,
            methods=['bonferroni', 'fdr_bh'],
            pvalcalc='fisher_scipy_stats'
            )
        results = goeaobj.run_study_nts(study_genes)

        sig_results = [r for r in results if r.p_uncorrected < 0.05] # TODO: May want to parse by minimum gene count as well
        print(len(sig_results))
        for i, r in enumerate(sig_results):
            if i > 10:
                break
            info = [str(s) for s in [r.GO, r.goterm.name, round(r.p_uncorrected, 5), round(r.p_fdr_bh, 5), r.ratio_in_study[0], r.ratio_in_study[1]]]
            print(info)
            print()

    @staticmethod
    def get_background_genes(background_genes_path: Path) -> None:
        # NOTE: This is a WIP
        # Currenty instructions are:
        # Go to: https://www.ncbi.nlm.nih.gov/gene
        # For Rousettus search: "9407"[Taxonomy ID] AND alive[property] AND genetype protein coding[Properties]
        # Send to file (must be tabular)
        pass
    
    @staticmethod
    def extract_background_genes(background_genes_path: Path) -> set[str]:
        background_genes = set()
        with background_genes_path.open() as inhandle:
            reader_iterator = reader(inhandle, delimiter="\t")
            header = next(reader_iterator)
            for line in reader_iterator:
                gene_id = int(line[2])
                background_genes.add(gene_id)
        return background_genes

    @staticmethod
    def set_symbol_id_mapper(background_genes_path: Path) -> dict[int]:
        symbol_id_mapper = {}
        with background_genes_path.open() as inhandle:
            reader_iterator = reader(inhandle, delimiter="\t")
            header = next(reader_iterator)
            for line in reader_iterator:
                gene_id = int(line[2])
                symbol = line[5]
                symbol_id_mapper[symbol] = gene_id
        return symbol_id_mapper
    
    @staticmethod
    def extract_study_genes(study_genes_path: Path) -> set[str]:
        study_genes = set()
        with study_genes_path.open() as inhandle:
            reader_iterator = reader(inhandle, delimiter=",")
            header = next(reader_iterator)
            for line in reader_iterator:
                gene = line[0]
                group = int(line[1])
                # if group != 12:
                #     continue
                study_genes.add(gene)
        return study_genes
    
    @staticmethod
    def convert_study_gene_symbols(study_genes: set[str], symbol_id_mapper: dict[str]) -> set[str]:
        converted_study_genes = set()
        for gene in study_genes:
            try:
                converted_study_genes.add(symbol_id_mapper[gene])
            except KeyError:
                continue
        return converted_study_genes

    @staticmethod
    def check_gene2go_id_count(gene2go_path: Path, taxid: str) -> None:
        gene_ids = set()
        with gene2go_path.open() as inhandle:
            reader_iterator = reader(inhandle, delimiter="\t")
            header = next(reader_iterator)
            for line in reader_iterator:
                if line[0] != taxid:
                    continue
                gene_id = line[1]
                gene_ids.add(gene_id)
        print(len(gene_ids))

if __name__ == "__main__":
    gm = GoatoolsManager()
    gm.setup()
    gm.run()