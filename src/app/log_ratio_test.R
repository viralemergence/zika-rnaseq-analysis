suppressMessages(library(DESeq2))

args = commandArgs(trailingOnly=TRUE)

gene_count_path <- args[1]
sample_metadata_path <- args[2]
lrt_results_outpath <- args[3]

# Load data
gene_counts=t(read.csv(gene_count_path, sep=",", head=T, row.names="X"))
sample_metadata=read.csv(sample_metadata_path, sep=",", head=T, row.names="Sample.ID")

# Define metadata factors for modelling
sample_metadata$Time <- factor(sample_metadata$Time)
sample_metadata$Virus <- factor(sample_metadata$Virus)

# Defining models
full_model <- ~ Virus + Time + Virus:Time
reduced_model <- ~ Virus + Time

# Creating DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = gene_counts,
                              colData = sample_metadata, 
                              design = full_model)

# Performing Likelihood Ratio Test to account for effect of time
dds_lrt_time <- DESeq(dds, test="LRT", reduced=reduced_model)

results_lrt <- results(dds_lrt_time)

# Write results to csv
write.csv(results_lrt, lrt_results_outpath, row.names=TRUE, quote=FALSE)