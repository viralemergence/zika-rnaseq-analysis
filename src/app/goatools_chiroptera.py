from collections import defaultdict
from csv import reader, DictWriter
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
        
        grouped_results = defaultdict(list)
        for r in sig_results:
            go_genes = [id_symbol_mapper[gene] for gene in r.study_items]
            if len(go_genes) < 2: # NOTE: May want to modulate threshold
                continue

            group_go_term = self.calculate_group_go_term(r.GO, godag, goslim_dag)
            
            info = {"go_id": r.GO,
                    "go_name": f'"{r.goterm.name}"',
                    "p-val": round(r.p_uncorrected, 5),
                    "genes": go_genes}
            grouped_results[group_go_term].append(info)

        grouped_results_highlights = list()
        for group_go_id, results in grouped_results.items():
            group_info = {"group_go_id": group_go_id,
                          "group_go_name": f'"{godag[group_go_id].name}"'}

            most_significant_result = sorted(results, key=lambda result: result["p-val"])[0]
            msr = most_significant_result
            most_significant_result_info = {"most_significant_go_id": msr["go_id"],
                                            "most_significant_go_name": msr["go_name"],
                                            "most_significant_p_val": msr["p-val"],
                                            "most_significant_study_genes": "|".join(msr["genes"])}
            group_info.update(most_significant_result_info)

            largest_result = sorted(results, key=lambda result: len(result["genes"]), reverse=True)[0]
            lr = largest_result
            largest_result_info = {"largest_go_id": lr["go_id"],
                                   "largest_go_name": lr["go_name"],
                                   "largest_p_val": lr["p-val"],
                                   "largest_study_genes": "|".join(lr["genes"])}
            group_info.update(largest_result_info)
            
            grouped_results_highlights.append(group_info)
        
        grouped_results_highlights = sorted(grouped_results_highlights, key=lambda result: result["most_significant_p_val"])
        
        with self.outpath.open("w") as outhandle:
            header = grouped_results_highlights[0].keys()
            writer = DictWriter(outhandle, fieldnames=header)
            writer.writeheader()
            writer.writerows(grouped_results_highlights)
    
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
                if group != 12:
                    continue
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

    @classmethod
    def calculate_group_go_term(cls, go_term: str, godag: GODag, goslim_dag: GODag) -> str:
        direct_ancestors, all_ancestors = mapslim(go_term, godag, goslim_dag)
        if direct_ancestors != all_ancestors:
            pass
            # print(f"{go_term} has different direct vs all ancestors")
        if len(all_ancestors) == 0:
            return "NONE"
        return cls.get_lowest_level_ancestor(godag, all_ancestors)

    @staticmethod
    def get_lowest_level_ancestor(godag: GODag, ancestors: set[str]) -> str:
        ancestor_levels = dict()
        for ancestor in ancestors:
            ancestor_levels[ancestor] = godag[ancestor].level
        return max(ancestor_levels, key= lambda x: ancestor_levels[x])

if __name__ == "__main__":
    taxon_id = 9407
    gm = GoatoolsManager(taxon_id)
    gm.setup()
    gm.run()