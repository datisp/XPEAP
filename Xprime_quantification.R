print("loading required R packages...")
library(tidyverse)
library(DESeq2)
library(RColorBrewer)
library(gplots)
library(Biostrings)

args <- commandArgs()
print(args)

strain1_rep1 <- args[6]
strain1_rep2 <- args[7]
strain1_rep3 <- args[8]
strain2_rep1 <- args[9]
strain2_rep2 <- args[10]
strain2_rep3 <- args[11]
cutoff_counts <- as.numeric(args[12])
path_fasta <- args[13]

### import bed files with read counts for Xprime ends
print("importing Xprime bed files...")
temp_list <- list.files(path = "./Xprime_analysis/wig2bed/coverage_Xprime", pattern = "*.bed", 
               full.names = TRUE)

read_id <- function(flnm) {
    read.csv(flnm, header = FALSE, sep = "\t") %>% 
        mutate(filename = flnm)
}

bed_imported <-
    list.files(path = "./Xprime_analysis/wig2bed/coverage_Xprime", pattern = "*.bed", 
               full.names = TRUE) %>% 
    map_df(~read_id(.))


### positive numbers only
print("converting all counts to positive numbers...")
bed_imported <- bed_imported %>%
  mutate(
    V5 = sqrt(V5 ^2)
  )

bed_imported <- bed_imported %>%
  filter(V5 >= cutoff_counts)

### define variables
seq_genome <- readDNAStringSet(path_fasta)
chrom_name <- NULL
chrom_length <- 0
for (i in seq_along(names(seq_genome))) {
  chrom_name[i] <- strsplit(names(seq_genome[i])," ")[[1]][1]
  chrom_length[i] <- width(seq_genome[i])
}

chrom_name
chrom_length


### create empty tibble for each forward counts and reverse counts. contains all chromosomes with corresponding position of nucleotides
print("generating nucleotide wise template files...")
# forward strands
gene_quanti_forward <- tibble(
  chromosome = character(),
  strand = character(),
  position = double()
)

for (i in seq_along(chrom_name)){
  tbl_tmp <- tibble(
    chromosome = chrom_name[i],
    strand = "+",
    position = 1:chrom_length[i]
  )
  gene_quanti_forward <- rbind(gene_quanti_forward, tbl_tmp)
}

# reverse strands
gene_quanti_reverse <- tibble(
  chromosome = character(),
  strand = character(),
  position = double()
)

for (i in seq_along(chrom_name)){
  tbl_tmp <- tibble(
    chromosome = chrom_name[i],
    strand = "-",
    position = 1:chrom_length[i]
  )
  gene_quanti_reverse <- rbind(gene_quanti_reverse, tbl_tmp)
}


# create lists containing all bed file names from + / - strand
print("generating file list for each strain...")
files_for <- c(
	paste0("./Xprime_analysis/wig2bed/coverage_Xprime/",strain1_rep1,"_trimmed_forward.wig.bed"),
	paste0("./Xprime_analysis/wig2bed/coverage_Xprime/",strain1_rep2,"_trimmed_forward.wig.bed"),
	paste0("./Xprime_analysis/wig2bed/coverage_Xprime/",strain1_rep3,"_trimmed_forward.wig.bed"),
	paste0("./Xprime_analysis/wig2bed/coverage_Xprime/",strain2_rep1,"_trimmed_forward.wig.bed"),
	paste0("./Xprime_analysis/wig2bed/coverage_Xprime/",strain2_rep2,"_trimmed_forward.wig.bed"),
	paste0("./Xprime_analysis/wig2bed/coverage_Xprime/",strain2_rep3,"_trimmed_forward.wig.bed")
)

files_rev <- c(
        paste0("./Xprime_analysis/wig2bed/coverage_Xprime/",strain1_rep1,"_trimmed_reverse.wig.bed"),
        paste0("./Xprime_analysis/wig2bed/coverage_Xprime/",strain1_rep2,"_trimmed_reverse.wig.bed"),
        paste0("./Xprime_analysis/wig2bed/coverage_Xprime/",strain1_rep3,"_trimmed_reverse.wig.bed"),
        paste0("./Xprime_analysis/wig2bed/coverage_Xprime/",strain2_rep1,"_trimmed_reverse.wig.bed"),
        paste0("./Xprime_analysis/wig2bed/coverage_Xprime/",strain2_rep2,"_trimmed_reverse.wig.bed"),
        paste0("./Xprime_analysis/wig2bed/coverage_Xprime/",strain2_rep3,"_trimmed_reverse.wig.bed")
)

files_for
files_rev

print("left joining all counts of all replicates and chromosomes/plasmids...")
# create empty tibbles for following loops, one for +/- each
prime_quanti_for <- tibble(
  chromosome = character(),
  strand = character(),
  position = double(),
  V5.x = double(),
  V5.y = double(),
  V5.x.x = double(),
  V5.y.y = double(),
  V5.x.x.x = double(),
  V5.y.y.y = double()
)

prime_quanti_rev <- tibble(
  chromosome = character(),
  strand = character(),
  position = double(),
  V5.x = double(),
  V5.y = double(),
  V5.x.x = double(),
  V5.y.y = double(),
  V5.x.x.x = double(),
  V5.y.y.y = double()
)

# join all Xprime counts of all bed - files with previously created template tibble which contains all nucleotide positions.
# one for + / - strand
# replicates and strains columnwise

for(c in seq_along(chrom_name)) {
tmp_bed <- filter(gene_quanti_forward, chromosome == chrom_name[c])
  for (i in seq_along(files_for)) {
      tmp <- filter(bed_imported, filename == files_for[i], V1 == chrom_name[c]) %>%
        select(V2, V5)
      tmp_bed <- left_join(tmp_bed, tmp, by = c("position" = "V2")) 
  }
prime_quanti_for <- rbind(prime_quanti_for, tmp_bed)
}

for(c in seq_along(chrom_name)) {
tmp_bed <- filter(gene_quanti_reverse, chromosome == chrom_name[c])
  for (i in seq_along(files_rev)) {
      tmp <- filter(bed_imported, filename == files_rev[i], V1 == chrom_name[c]) %>%
        select(V2, V5)
      tmp_bed <- left_join(tmp_bed, tmp, by = c("position" = "V2")) 
  }
prime_quanti_rev <- rbind(prime_quanti_rev, tmp_bed)
}

# bind prime_quanti_for and _rev to one tibble
prime_quanti_full <- rbind(prime_quanti_for, prime_quanti_rev)

# discard all NA-rows and substitute all remaining NA with 1 (pseudo count) to enable downstream DESeq2 analysis
prime_quanti_full <- filter(prime_quanti_full, !((is.na(V5.x)) & (is.na(V5.y)) & (is.na(V5.x.x)) & (is.na(V5.y.y)) & (is.na(V5.x.x.x)) & (is.na(V5.y.y.y)))  ) %>%
  replace_na(list(V5.x = 0, V5.y = 0, V5.x.x = 0, V5.y.y = 0, V5.x.x.x = 0, V5.y.y.y = 0))


quantification_full <- prime_quanti_full %>%
	dplyr::rename(
		"mutant_biorep1" = "V5.x",
		"mutant_biorep2" = "V5.y",
		"mutant_biorep3" = "V5.x.x",
		"wt_biorep1" = "V5.y.y",
		"wt_biorep2" = "V5.x.x.x",
		"wt_biorep3" = "V5.y.y.y"
	)

head(quantification_full)

print("exporting Xprime quantification file to ./Xprime_analysis/Xprime_DESeq/quantification_full.csv...")
write_tsv(quantification_full, "./Xprime_analysis/Xprime_DESeq/quantification_full.csv")
