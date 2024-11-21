from csv import reader
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.base import download_go_basic_obo, dnld_file
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from goatools.obo_parser import GODag
from pathlib import Path

class GoatoolsManager:
    def __init__(self) -> None:
        self.background_genes_path = Path("/src/data/pydeseq2/goatools/ncbi_gene_results.txt")
        self.obo_path = Path("/src/data/pydeseq2/goatools/go-basic.obo")
        self.gene2go_ftp_path = "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz"
        self.gene2go_path = Path("/src/data/pydeseq2/goatools/gene2go.txt")
        self.gene_list = {"PDLIM1", "SNAPC1", "ISG15", "PAQR5", "LOC107513409",
                          "MTSS1", "TNFAIP6", "LOC118611464", "IFI6", "SEPTIN3",
                          "LOC107512537", "FAM110A", "TUBD1", "MX1", "LOC118606763",
                          "LOC118612456", "ID3", "SPSB1"} #NOTE: Temp for testing

    def setup(self) -> None:
        download_go_basic_obo(str(self.obo_path))
        dnld_file(self.gene2go_ftp_path, str(self.gene2go_path))
        self.get_background_genes(self.background_genes_path)
        
    def run(self) -> None:
        background_genes = self.extract_background_genes(self.background_genes_path)
        symbol_id_mapper = self.set_symbol_id_mapper(self.background_genes_path)
        study_genes = self.convert_study_gene_symbols(self.gene_list, symbol_id_mapper)

        godag = GODag(self.obo_path) # NOTE: Check this

        objanno = Gene2GoReader(str(self.gene2go_path), taxids=[9407])
        ns2assoc = objanno.get_ns2assc()
        goeaobj = GOEnrichmentStudyNS(
            background_genes,
            ns2assoc,
            godag,
            propagate_counts=False,
            alpha=0.5,
            methods=['fdr_bh']
            ) # NOTE: DO NOT LEAVE AT CURRENT ALPHA VALUE!!!
        
        results = goeaobj.run_study(study_genes)
        
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
    def convert_study_gene_symbols(study_genes: set[str], symbol_id_mapper: dict[str]) -> set[str]:
        converted_study_genes = set()
        for gene in study_genes:
            try:
                converted_study_genes.add(symbol_id_mapper[gene])
            except KeyError:
                continue
        return converted_study_genes

if __name__ == "__main__":
    gm = GoatoolsManager()
    gm.setup()
    gm.run()