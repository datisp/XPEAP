#!/bin/bash

echo "converting BED to sorted BED"
filename_input="deseq_"$strain1"_vs_"$strain2"_Xprime_ends.bed"
filename_sorted="deseq_"$strain1"_vs_"$strain2"_Xprime_ends_sorted.bed"
filename_merged="deseq_"$strain1"_vs_"$strain2"_Xprime_ends_sorted_merged.bed"
filename_intersection=$strain1"_vs_"$strain2"_intersection_ggf_Xprime_ends.bed"

bedtools sort -i ./Xprime_analysis/Xprime_DESeq/$filename_input > ./Xprime_analysis/Xprime_DESeq/$filename_sorted

echo "merging all positions within a distance of 3 nucleotides"
echo "mean value of score (log2FC) ist computed and saved in colum 5"
echo "strandedness is forced"
echo "sum of all baseMeans is computed and saved in colum 6"
bedtools merge -s -d 3 -c 5,7 -o mean,sum -i ./Xprime_analysis/Xprime_DESeq/$filename_sorted > ./Xprime_analysis/Xprime_DESeq/$filename_merged

echo "performing intersection analysis between annotation.gff and deseq_3prime_merged.bed"
bedtools intersect -wo -a ./Xprime_analysis/annotation.bed -b ./Xprime_analysis/Xprime_DESeq/$filename_merged > ./Xprime_analysis/$filename_intersection
echo "finished"
