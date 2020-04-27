#!/bin/bash

### READemption pipeline: X'end anaylsis

echo "reademption_complete.sh started..."

# activate READemption conda env
source activate $conda_reademption
#bash
#source activate rnaseq

##reademption create $dir

date > ./$dir/log_align
date > ./$dir/runtime

echo "copying trimmed read files..."
##cp ./reads_trimmed/*trimmed.fq ./$dir/input/reads/
##cp $path_genomic_fasta ./$dir/input/reference_sequences/ref_seq.fa
##cp $path_gff ./$dir/input/annotations/annotation.gff

echo "running READemption align..."
# READemption align
##reademption align -p $n_cores --progress --fastq ./$dir &>> ./$dir/log

echo "running READemption gene quantification"
# READemption gene quanti
##date > ./$dir/log_gene_quanti
##reademption gene_quanti -p $n_cores --features CDS,tRNA,rRNA,tmRNA,ncRNA,sRNA ./$dir &>> ./$dir/log_gene_quanti


echo "DESeq2..."
##date > ./$dir/log_deseq
##reademption deseq --libs $input_libs --conditions $input_conditions ./$dir &>> ./$dir/log_deseq

# READemption visualisation
##reademption viz_deseq ./$dir
##reademption coverage -p $n_cores ./$dir

# rename coverage folder
##mv ./$dir/output/coverage ./$dir/output/coverage_full

# generate new folder structure for Xprime mapping
##mkdir -p ./$dir/output/coverage/{coverage-raw,coverage-tnoar_mil_normalized,coverage-tnoar_min_normalized}

echo "running X'mapping for replicates..."
# generate Xprime mapping
##reademption coverage -p $n_cores --coverage_style $coverage_style ./$dir
##mv ./$dir/output/coverage ./$dir/output/coverage_Xprime

#remove dublicated read files to save disc space
##rm ./$dir/input/reads/*
###echo "READemption analysis for replicate files finished"

#================================================================================#
#================================================================================#

echo "merging trimmed read files..."
# folder structure merged files
dir_merged=$dir"_merged"
filetype=".fq"
reademption create $dir_merged
cp ./reads_trimmed/*trimmed$filetype ./$dir_merged/input/reads/

# merge read files
dir_reads=$dir_merged"/input/reads"
cat {./$dir_reads/$strain1_rep1"_trimmed"$filetype,./$dir_reads/$strain1_rep2"_trimmed"$filetype,./$dir_reads/$strain1_rep3"_trimmed"$filetype} > ./$dir_reads/$strain1"_merged.fq"
cat {./$dir_reads/$strain2_rep1"_trimmed"$filetype,./$dir_reads/$strain2_rep2"_trimmed"$filetype,./$dir_reads/$strain2_rep3"_trimmed"$filetype} > ./$dir_reads/$strain2"_merged.fq"
rm -r ./$dir_reads/*trimmed*
cp $path_genomic_fasta ./$dir_merged/input/reference_sequences/ref_seq.fa
cp $path_gff ./$dir_merged/input/annotations/annotation.gff

echo "running READemption align..."
# READemption align
reademption align -p $n_cores --progress --fastq ./$dir_merged &>> ./$dir_merged/log

echo "running READemption coverage..."
#READemption converage
reademption coverage -p $n_cores ./$dir_merged

# rename coverage folder
mv ./$dir_merged/output/coverage ./$dir_merged/output/coverage_full_merged

# generate new folder structure for Xprime mapping
mkdir -p ./$dir_merged/output/coverage/{coverage-raw,coverage-tnoar_mil_normalized,coverage-tnoar_min_normalized}

# generate 3prime mapping
reademption coverage -p $n_cores --coverage_style $coverage_style ./$dir_merged
mv ./$dir_merged/output/coverage ./$dir_merged/output/coverage_Xprime_merged

#remove dublicated read files to save disc space
rm ./$dir_merged/input/reads/*

echo "reademption_complete.sh finished."
date >> ./$dir/runtime
