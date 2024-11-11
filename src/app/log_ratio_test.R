print("Starting")
library(DESeq2)

gene_count_path <- "/home/alexander.brown/zika-rnaseq-analysis/src/data/pydeseq2/gene_counts_R_input.csv"
sample_metadata_path <- "/home/alexander.brown/zika-rnaseq-analysis/src/data/pydeseq2/sample_metadata_R_input.csv"

gene_counts=read.csv(gene_count_path, sep=",", head=T, row.names="X")
sample_metadata=read.csv(sample_metadata_path, sep=",", head=T, row.names="Sample.ID")

sample_metadata$Time <- factor(sample_metadata$Time)
sample_metadata$Virus <- factor(sample_metadata$Virus)
#head(sample_metadata, n=c(20,6))

full_model <- ~ Virus + Time + Virus:Time
reduced_model <- ~ Virus + Time

print(full_model)

print("Done")