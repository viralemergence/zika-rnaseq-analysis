library(DEGreport)

args = commandArgs(trailingOnly=TRUE)

gene_count_path <- args[1]
sample_metadata_path <- args[2]

# Load data
gene_counts=t(read.csv(gene_count_path, sep=",", head=T, row.names="X"))
sample_metadata=read.csv(sample_metadata_path, sep=",", head=T, row.names="Sample.ID")

pdf(file="/src/data/pydeseq2/test.pdf")
clusters <- degPatterns(gene_counts, metadata=sample_metadata, time="Time", col="Virus", plot=T)
dev.off()