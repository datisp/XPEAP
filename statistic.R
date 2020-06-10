library(tidyverse)

args <- commandArgs()
print(args)

merged_bed <- read_tsv(paste0("./Xprime_analysis/Xprime_DESeq/deseq_",args[6],"_vs_",args[7],"_Xprime_ends_sorted_merged.bed"), col_names = F) %>%
  dplyr::rename("chromosome" = "X1",
                "start" = "X2",
                "end" = "X3",
                "strand" = "X4",
                "log2FoldChange" = "X5",
		"baseMean" = "X6") %>%
  mutate(length = end - start)

# group all detected and valid X'-ends in clusters according to their proximity
cluster_var <- 1
windowsize <- 25

merged_bed_cluster <- merged_bed %>%
  add_column(
    cluster_ID = 1,
    significance = NA,
    Xprime_end_ID = 1:nrow(merged_bed)
  )

for (i in (seq_along(1:(nrow(merged_bed_cluster)-1)))) {
  if ((merged_bed_cluster$end[i] + windowsize) >= merged_bed_cluster$start[i+1]) {
    merged_bed_cluster$cluster_ID[i+1] <- cluster_var
  } else {
    cluster_var <- (cluster_var + 1)
    merged_bed_cluster$cluster_ID[i+1] <- cluster_var
  }
}

# find maximum baseMean value in each cluster
# mark corresponding Xprime end as "major" Xprime end
cluster_no <- length(unique(merged_bed_cluster$cluster_ID))
for (i in seq_along(1:cluster_no)) {
  cluster_tmp <- filter(merged_bed_cluster, cluster_ID == i)
  coverage_max <- max(cluster_tmp$baseMean)
  cluster_max <- filter(cluster_tmp, baseMean == coverage_max)
  for (j in seq_along(1:nrow(merged_bed_cluster))) {
    if (merged_bed_cluster$Xprime_end_ID[j] == cluster_max$Xprime_end_ID) {
      merged_bed_cluster$significance[j] <- "major"
    } else {
      NULL
    }
  }
}

# mark all other Xprime ends as "minor"
for (i in (seq_along(1:nrow(merged_bed_cluster)))) {
  if (is.na(merged_bed_cluster$significance[i])) {
    merged_bed_cluster$significance[i] <- "minor"
  } else {
    NULL
  }
}

# filter all major Xprime ends and save output
merged_bed_major <- merged_bed_cluster %>%
	filter(significance == "major")
head(merged_bed_major)
print(paste0("writing BED file containing all MAJOR X'-ends to ./Xprime_analysis/Xprime_DESeq/deseq_",args[6],"_vs_",args[7],"_Xprime_ends_sorted_merged_major.bed"))
write_tsv(merged_bed_major, paste0("./Xprime_analysis/Xprime_DESeq/deseq_",args[6],"_vs_",args[7],"_Xprime_ends_sorted_merged_major.bed"), col_names = F)

# basic summary plots
p_histogram_3ends_differential_log2FC <- merged_bed %>%
  ggplot(aes(log2FoldChange)) +
  geom_histogram() +
  labs(subtitle = "differential X'ends: count of all log2FC", x = "log2FC mutant vs. wildtype") +
  theme_bw()
p_histogram_3ends_differential_log2FC

p_histogram_3ends_length <- merged_bed %>%
  ggplot(aes(length)) +
  geom_histogram(binwidth = 3) +
  labs(subtitle = "differential X'ends: count of all lengths", x = "length of X'end after bedtools:merge [nt]") +
  theme_bw()
p_histogram_3ends_length

p_histogram_3ends_length_zoom <- merged_bed %>%
  filter(length <= 50) %>%
  ggplot(aes(length)) +
  geom_histogram(binwidth = 3) +
  labs(subtitle = "differential X'ends: count of all lengths - zoomed in", x = "length of X'end after bedtools:merge [nt]") +
  theme_bw()
p_histogram_3ends_length_zoom

# read intersection file
intersections <- read_tsv(paste0("./Xprime_analysis/",args[6],"_vs_",args[7],"_intersection_ggf_Xprime_ends.bed"), col_names = FALSE) %>%
  dplyr::rename("chromosome" = "X1",
                "feature_start" = "X2",
                "feature_end" = "X3",
                "attribute_name" = "X4",
                "feature" = "X5",
                "start" = "X7",
                "end" = "X8",
                "strand" = "X9",
                "log2FoldChange" = "X10") %>%
  select(chromosome, feature, start, end, log2FoldChange, feature_start, feature_end, attribute_name)

# define UTR overlap: all X'ends that do not intersect with provided gff have to be originated in UTRs
UTR_overlap <- left_join(merged_bed, intersections, by = "start") %>%
  filter(is.na(feature)) %>%
  dplyr::rename(
    "chromosome" = "chromosome.x",
    "end" = "end.x",
    #"strand" = "strand.x",
    "log2FoldChange" = "log2FoldChange.x") %>%
  select(chromosome, feature, start, end, log2FoldChange, feature_start, feature_end) %>%
  mutate(feature = "UTR") %>%
  add_column(attribute_name = "unknown")

# rbind intersections and UTRs
intersections_full <- rbind(intersections, UTR_overlap) %>%
  add_column(change = "unknown") %>%
  group_by(attribute_name) %>%
  mutate(
    count = n(),
    feature_length = feature_end - feature_start,
    count_per_kb = (count/feature_length)*1000)

for (i in seq_along(1:nrow(intersections_full))) {
      if (intersections_full$log2FoldChange[i] > 0) {
     intersections_full$change[i] <- "mutant enriched X'ends"
   } else {
     intersections_full$change[i] <- "wildtype enriched X'ends"
   }
}

# group features: sRNA + tmRNA -> ncRNA
for (i in seq_along(1:nrow(intersections_full))) {
  if (intersections_full$feature[i] == "sRNA") {
    intersections_full$feature[i] <- "ncRNA"
  } else if (intersections_full$feature[i] == "tmRNA") {
    intersections_full$feature[i] <- "ncRNA"
  } else {
    NULL
  }
}

write_csv2(intersections_full, "./Xprime_analysis/intersections_full.csv")
 
#### plots 
alpha_boxplot <- 0

p_histogram_3ends_differential <-  intersections_full %>%
  ggplot(aes(feature, fill = change)) +
  geom_histogram(aes(), stat = "count", position = position_dodge(), color = "black") +
  theme_bw()
p_histogram_3ends_differential

### filter intersection list genewise
intersections_genewise <- intersections_full %>%
  group_by(change) %>%
  distinct(attribute_name, feature, count, feature_length, count_per_kb) 

### plots
p_boxplot_3ends_abs <- ggplot(data = filter(intersections_genewise, feature != "UTR"), aes(feature, count)) +
  geom_jitter(aes(color = feature), alpha = 0.5, width = 0.3) +
  geom_boxplot(outlier.alpha = FALSE, alpha = alpha_boxplot) +
  labs(title = "absolute counts", y = "number of differential X'ends per gene") +
  theme_bw()
p_boxplot_3ends_abs

p_boxplot_3ends_rel <- ggplot(data = filter(intersections_genewise, feature != "UTR"), aes(feature, count_per_kb)) +
  geom_jitter(aes(color = feature), alpha = 0.5, width = 0.3) +
  geom_boxplot(outlier.alpha = FALSE, alpha = alpha_boxplot) +
  labs(title = "counts per kb", y = "number of differential X'ends per gene per kb") +
  theme_bw()
p_boxplot_3ends_rel

p_boxplot_3ends_log2FC <- ggplot(intersections_full, aes(feature, log2FoldChange)) +
  geom_jitter(aes(color = feature), alpha = 0.5, width = 0.3) +
  geom_boxplot(outlier.alpha = FALSE, alpha = alpha_boxplot) +
  labs(title = "all differential X'ends") +
  theme_bw()
p_boxplot_3ends_log2FC

p_histogram_3ends_differential_featurewise_abs <- ggplot(intersections_full, aes(log2FoldChange)) +
  geom_histogram(aes(fill = feature), color = "black") +
  scale_x_continuous(breaks = seq(-10,10,5)) +
  facet_wrap(~ feature, nrow = 2) +
  labs(title = "mutant vs. wt, detected differential X'ends (absolute counts)") +
  theme_bw()

p_histogram_3ends_differential_featurewise_rel <- ggplot(intersections_full, aes(log2FoldChange)) +
  geom_histogram(aes(y = ..density.. , fill = feature), color = "black") +
  scale_x_continuous(breaks = seq(-10,10,5)) +
  facet_wrap(~ feature, nrow = 2) +
  labs(title = "mutant vs. wt, detected differential X'ends (relative counts)") +
  theme_bw()

pdf("./Xprime_analysis/figures_statistic.pdf", width = 6, height = 4)
p_histogram_3ends_differential_log2FC
p_histogram_3ends_length
p_histogram_3ends_length_zoom
p_histogram_3ends_differential
p_boxplot_3ends_abs
p_boxplot_3ends_rel
p_histogram_3ends_differential_featurewise_abs
p_histogram_3ends_differential_featurewise_rel
dev.off()

svg("./Xprime_analysis/figures_statistic%02d.svg", width = 6, height = 4)
p_histogram_3ends_differential_log2FC
p_histogram_3ends_length
p_histogram_3ends_length_zoom
p_histogram_3ends_differential
p_boxplot_3ends_abs
p_boxplot_3ends_rel
p_histogram_3ends_differential_featurewise_abs
p_histogram_3ends_differential_featurewise_rel
dev.off()


print("figures written to ./Xprime_analysis/")


#================================================================================#
stats_general <- tibble("category" = 1:6,
                    "number" = 0)

### all differential X'ends (merged nucleotides)
stats_general$category[1] <- "all differential X'ends"
stats_general$number[1] <- nrow(merged_bed)
stats_general$category[2] <- "mutant enriched X'ends"
stats_general$number[2] <- merged_bed %>%
  filter(log2FoldChange > 0) %>%
  nrow()
stats_general$category[3] <- "wildtype enriched X'ends"
stats_general$number[3] <- merged_bed %>%
  filter(log2FoldChange <= 0) %>%
  nrow()


### intersections: bedtools intersect counts
# number of X'ends per feature
stats_general$category[4] <- "all differential X'ends after bedtools:intersect analysis"
stats_general$number[4] <- nrow(intersections_full)

stats_intersection_featurewise <- intersections_full %>%
  group_by(change, feature) %>%
  summarize(
    count = n()
  )

# number of genes with differential X'ends. since there is no valid annotatioin for separate UTRs, all UTRs count as ONE feature!
stats_differential_genes <- intersections_full %>%
  group_by(change, attribute_name) %>%
  summarize(
    count = n()
  ) %>%
  nrow()

stats_general$category[5] <- "genes with mutant enriched X'ends (all UTRs count as 1)"
stats_general$number[5] <- intersections_full %>%
  filter(change == "mutant enriched X'ends") %>%
  group_by(attribute_name) %>%
  summarize(
    count = n()
  ) %>%
  nrow()

stats_general$category[6] <- "genes with wildtype enriched X'ends (all UTRs count as 1)"
stats_general$number[6] <- intersections_full %>%
  filter(change == "wildtype enriched X'ends") %>%
  group_by(attribute_name) %>%
  summarize(
    count = n()
  ) %>%
  nrow()

write_csv2(stats_general, "./Xprime_analysis/stats_general.csv")
print("general statistic information written to ./Xprime_analysis/stats_general.csv")
