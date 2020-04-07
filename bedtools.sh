#!/bin/bash

echo "converting BED to sorted BED"
filename_input="deseq_"$condition1"_vs_"$condition2"_3prime_ends.bed"
filename_sorted="deseq_"$condition1"_vs_"$condition2"_3prime_ends_sorted.bed"
filename_merged="deseq_"$condition1"_vs_"$condition2"_3prime_ends_sorted_merged.bed"
filename_intersection=$condition1"_vs_"$condition2"_intersection_ggf_3prime_ends.bed"

bedtools sort -i ./3prime_analysis/3prime_DESeq/$filename_input > ./3prime_analysis/3prime_DESeq/$filename_sorted

echo "merging all positions within a distance of 3 nucleotides"
echo "mean value of score (log2FC) ist computed and saved in colum 5"
echo "strandedness is forced"
bedtools merge -s -d 3 -c 5 -o mean -i ./3prime_analysis/3prime_DESeq/$filename_sorted > ./3prime_analysis/3prime_DESeq/$filename_merged

echo "performing intersection analysis between annotation.gff and deseq_3prime_merged.bed"
bedtools intersect -wo -a ./3prime_analysis/annotation.bed -b ./3prime_analysis/3prime_DESeq/$filename_merged > ./3prime_analysis/$filename_intersection
echo "finished"
