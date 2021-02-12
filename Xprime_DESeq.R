library(tidyverse)
library(DESeq2)
library(RColorBrewer)
library(gplots)

args <- commandArgs()
print(args)

### DESeq2
print("starting DESeq2 analysis...")
count_table_raw <- read_tsv("./Xprime_analysis/Xprime_DESeq/quantification_full.csv")
position_metadata <- count_table_raw[3]
count_table <- round(count_table_raw[,4:9])
libraries <- c(args[6], args[7], args[8], args[9], args[10], args[11])
conditions <- c(args[12],args[12],args[12],args[13],args[13],args[13])
colnames(count_table) <- libraries
samples <- data.frame(row.names=libraries, condition=conditions, lib=libraries)
dds <- DESeqDataSetFromMatrix(countData=count_table, colData=samples, design=~condition)
dds <- DESeq(dds)

### plots
pdf(paste0("./Xprime_analysis/Xprime_DESeq/PCA_Xprime_ends_",args[12],"_vs_",args[13],".pdf"))
rld <- rlog(dds)
print(plotPCA(rld, intgroup=c('condition')))

svg(paste0("./Xprime_analysis/Xprime_DESeq/PCA_Xprime_ends_",args[12],"_vs_",args[13],"%02d.svg"))
rld <- rlog(dds)
print(plotPCA(rld, intgroup=c('condition')))

distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- with(colData(dds), paste(lib, sep=' : '))
hmcol <- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
heatmap.2(mat, trace='none', col = rev(hmcol), margin=c(13, 13))

dds_result <- results(dds, contrast=c('condition',args[12], args[13]))
prime_deseq <- cbind(count_table_raw, dds_result)
prime_deseq <- as.tibble(prime_deseq)

print(paste0("writing raw DESeq2 analysis file to ./Xprime_analysis/Xprime_DESeq/deseq_",args[12],"_vs_",args[13],"_Xprime_ends.csv ... "))
write_tsv(prime_deseq, paste0("./Xprime_analysis/Xprime_DESeq/deseq_",args[12],"_vs_",args[13],"_Xprime_ends.csv"))
print("check")

### cutoff variables
cutoff_log2FC <- as.numeric(args[14])
cutoff_padj <- as.numeric(args[15])

prime_deseq_filtered <- prime_deseq %>%
  filter(padj <= cutoff_padj) %>%
  filter(sqrt(log2FoldChange^2) >= cutoff_log2FC)

print(paste0("writing filtered DESeq2 analysis file to ./Xprime_analysis/Xprime_DESeq/deseq_",args[12],"_vs_",args[13],"_Xprime_ends_filtered.csv ... "))
write_tsv(prime_deseq_filtered, paste0("./Xprime_analysis/Xprime_DESeq/deseq_",args[12],"_vs_",args[13],"_Xprime_ends_filtered.csv"))
print("check")

### add identifier
prime_deseq_filtered <- prime_deseq_filtered %>%
  add_column(ID = 1:nrow(prime_deseq_filtered))

### convert deseq file to bed file
prime_bed <- prime_deseq_filtered %>%
  mutate(position_end = position +1) %>%
  select(chromosome, position, position_end, ID, log2FoldChange, strand, baseMean)

print("converting filtered DESeq2 analysis file to BED format...")
write_tsv(prime_bed, paste0("./Xprime_analysis/Xprime_DESeq/deseq_",args[12],"_vs_",args[13],"_Xprime_ends.bed"), col_names = FALSE)
print("check")

print("generating optimized BED file instead of GFF for further intersection analysis")
filename_deseq <- paste0("./",args[16],"/output/deseq/deseq_with_annotations/deseq_comp_",args[12],"_vs_",args[13],"_with_annotation_and_countings.csv")
annotation_bed <- read_tsv(filename_deseq, skip = 2) %>%
  filter(`Orientation of counted reads relative to the strand location of the annotation` == "sense") %>%
  select("Sequence name", "Start", "End", "Attributes", "Feature")
write_tsv(annotation_bed, "./Xprime_analysis/annotation.bed", col_names = FALSE)
print("finished")
