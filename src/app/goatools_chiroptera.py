from collections import defaultdict
from csv import reader
from goatools.anno.genetogo_reader import Gene2GoReader # type: ignore
from goatools.base import dnld_file # type: ignore
from goatools.go_enrichment import GOEnrichmentStudy # type: ignore
from goatools.mapslim import mapslim # type: ignore
from goatools.obo_parser import GODag # type: ignore
from pathlib import Path

class GoatoolsManager:
    def __init__(self, taxon_id: int) -> None:
        # Background and study genes for enrichment analysis
        self.background_genes_path = Path("/src/data/go_analysis/ncbi_gene_results.txt")
        self.study_genes_path = Path("/src/data/pydeseq2/degpatterns/gene_clusters.csv")

        # Full and slim GO terms in obo format
        self.obo_ftp_path = "http://current.geneontology.org/ontology/go.obo"
        self.obo_path = Path("/src/data/go_analysis/go.obo")
        self.slim_obo_ftp_path = "https://current.geneontology.org/ontology/subsets/goslim_pir.obo"
        self.slim_obo_path = Path("/src/data/go_analysis/goslim-pir.obo")

        # gene2go and taxon specific gene2go info
        self.gene2go_ftp_path = "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz"
        self.gene2go_path = Path("/src/data/go_analysis/gene2go.txt")
        self.taxon_id = taxon_id
        self.gene2go_taxon_path = Path(f"/src/data/go_analysis/gene2go_{self.taxon_id}.txt")

        # Final results path
        self.outpath = Path("/src/data/go_analysis/goatools_results.csv")

    def setup(self) -> None:
        if not self.gene2go_taxon_path.is_file():
            dnld_file(self.gene2go_ftp_path, str(self.gene2go_path))
            self.filter_gene2go_by_taxon_id(self.gene2go_path, self.gene2go_taxon_path, self.taxon_id)
        dnld_file(self.obo_ftp_path, str(self.obo_path))
        dnld_file(self.slim_obo_ftp_path, str(self.slim_obo_path))
        self.get_background_genes(self.background_genes_path)
        
    @staticmethod
    def filter_gene2go_by_taxon_id(gene2go_path: Path, gene2go_taxon_path: Path, taxon_id: int) -> None:
        # Creating a taxon specific gene2go file dramatically increases Gene2GoReader load time
        with gene2go_path.open() as inhandle, gene2go_taxon_path.open("w") as outhandle:
            reader_iterator = reader(inhandle, delimiter="\t")
            header = next(reader_iterator)
            outhandle.write("\t".join(header) + "\n")
            for line in reader_iterator:
                tax_id = int(line[0])
                if tax_id != taxon_id:
                    continue
                outhandle.write("\t".join(line) + "\n")

    @staticmethod
    def get_background_genes(background_genes_path: Path) -> None:
        # NOTE: This is a WIP
        # Currenty instructions are:
        # Go to: https://www.ncbi.nlm.nih.gov/gene
        # For Rousettus search: "9407"[Taxonomy ID] AND alive[property] AND genetype protein coding[Properties]
        # Send to file (must be tabular)
        pass

    def run(self) -> None:
        background_genes = self.extract_background_genes(self.background_genes_path)
        symbol_id_mapper = self.set_symbol_id_mapper(self.background_genes_path)
        id_symbol_mapper = {id: symbol for symbol, id in symbol_id_mapper.items()}
        study_genes = self.extract_study_genes(self.study_genes_path)
        study_genes = self.convert_study_gene_symbols(study_genes, symbol_id_mapper)

        godag = GODag(str(self.obo_path), optional_attrs=["relationship"])
        goslim_dag = GODag(str(self.slim_obo_path), optional_attrs=["relationship"])

        objanno = Gene2GoReader(str(self.gene2go_taxon_path), taxids=[self.taxon_id])
        ns2assoc = objanno.get_ns2assc()
        id2gos = defaultdict(set)
        for _, associations in ns2assoc.items():
            for gene_id, go_terms in associations.items():
                id2gos[gene_id].update(go_terms)

        goeaobj = GOEnrichmentStudy(
            background_genes,
            id2gos,
            godag,
            propogate_counts=True,
            methods=['bonferroni', 'fdr_bh'],
            pvalcalc='fisher_scipy_stats'
            )
        results = goeaobj.run_study_nts(study_genes)

        sig_results = [r for r in results if r.p_uncorrected < 0.05]
        most_genes_per_term_results = sorted(sig_results, key=lambda r: r.ratio_in_study[0], reverse=True)
        sig_results = most_genes_per_term_results
        
        with self.outpath.open("w") as outhandle:
            for r in sig_results:
                go_genes = [id_symbol_mapper[gene] for gene in r.study_items]
                if len(go_genes) < 2: # NOTE: May want to modulate threshold
                    continue
                direct_ancestors, all_ancestors = mapslim(r.GO, godag, goslim_dag)
                if direct_ancestors != all_ancestors:
                    print(r.GO)
                if len(all_ancestors) == 0:
                    all_ancestors = ["NONE"]
                go_genes = "|".join(go_genes)
                all_ancestors = "|".join(all_ancestors)
                info = [r.GO, r.goterm.name, round(r.p_uncorrected, 5), round(r.p_fdr_bh, 5), r.ratio_in_study[0], r.ratio_in_study[1], go_genes, all_ancestors]
                info = ",".join([str(s) for s in info])
                outhandle.write(info + "\n")
    
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

if __name__ == "__main__":
    taxon_id = 9407
    gm = GoatoolsManager(taxon_id)
    gm.setup()
    gm.run()