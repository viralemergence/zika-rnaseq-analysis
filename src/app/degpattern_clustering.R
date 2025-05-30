library(DEGreport)

args = commandArgs(trailingOnly=TRUE)

gene_count_path <- args[1]
sample_metadata_path <- args[2]
degpatterns_results_outpath <- args[3]
degpatterns_figure_outpath <- args[4]

# Load data
gene_counts=t(read.csv(gene_count_path, sep=",", head=T, row.names="X"))
sample_metadata=read.csv(sample_metadata_path, sep=",", head=T, row.names="Sample.ID")

# Define metadata factors for clustering
sample_metadata$Time <- factor(sample_metadata$Time)
sample_metadata$Virus <- factor(sample_metadata$Virus)

pdf(file=degpatterns_figure_outpath)
clusters <- degPatterns(log2(gene_counts), metadata=sample_metadata, time="Time", col="Virus", scale=T, plot=T)
# Not sure if log2 is best here, but degpattern says input should be log2 normalized count matrix
dev.off()

# Write results to csv
write.csv(clusters$df, degpatterns_results_outpath, row.names=F, quote=F)