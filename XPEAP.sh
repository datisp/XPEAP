#!/bin/bash

echo "     ___          ___       ___          ___          ___   "
echo "    /__/|        /  /\     /  /\        /  /\        /  /\  "
echo "   |  |:|       /  /::\   /  /:/_      /  /::\      /  /::\ "
echo "   |  |:|      /  /:/\:\ /  /:/ /\    /  /:/\:\    /  /:/\:\ "
echo " __|__|:|     /  /:/~/://  /:/ /:/_  /  /:/~/::\  /  /:/~/:/"
echo "/__/::::\____/__/:/ /://__/:/ /:/ /\/__/:/ /:/\:\/__/:/ /:/ "
echo "   ~\~~\::::/\  \:\/:/ \  \:\/:/ /:/\  \:\/:/__\/\  \:\/:/  "
echo "    |~~|:|~~  \  \::/   \  \::/ /:/  \  \::/      \  \::/   "
echo "    |  |:|     \  \:\    \  \:\/:/    \  \:\       \  \:\   "
echo "    |  |:|      \  \:\    \  \::/      \  \:\       \  \:\  "
echo "    |__|/        \__\/     \__\/        \__\/        \__\/  "

# 2021-02-12 version 1.0.1

#================================================================================#

# please modify the content of the following variables:

# path to folder with raw data
export path_raw_data="/path/raw/data"

# path to gff and genomic.fasta
export path_gff="/path/gff/genomic.gff"
export path_genomic_fasta="/path/gff/genomic.fa"

export strain1="strain1"
export strain2="strain1"

# save filenames of each read file (without file type extension!)
export strain1_rep1="filename_strain1_rep1"
export strain1_rep2="filename_strain1_rep2"
export strain1_rep3="filename_strain1_rep3"
export strain2_rep1="filename_strain2_rep1"
export strain2_rep2="filename_strain2_rep2"
export strain2_rep3="filename_strain2_rep3"
export filetype=".fq"
export compression=".gz"

# path to fastqc bin file
export path_fastqc="/path/FastQC"

# names of conda environments
export conda_multiqc="name_conda_multiqc"
export conda_trim_galore="name_conda_trimming"
export conda_reademption="name_conda_rnaseq"
export conda_bedops="name_conda_bedops"

export cutoff_counts="10"
export cutoff_padj="0.05"
export cutoff_log2FC="1"
export cutoff_ratio="0.05"

# READemption settings:
export n_cores="24" # number of cores to be used for computation
export dir="name_dir" # name for READemption project folders
export coverage_style="last_base_only" # options are "first_base_only" and "last_base_only"

#================================================================================#

# READemption completed variables
export input_libs=$strain1_rep1"_trimmed"$filetype","$strain1_rep2"_trimmed"$filetype","$strain1_rep3"_trimmed"$filetype","$strain2_rep1"_trimmed"$filetype","$strain2_rep2"_trimmed"$filetype","$strain2_rep3"_trimmed"$filetype
export input_conditions=$strain1","$strain1","$strain1","$strain2","$strain2","$strain2


#================================================================================#
# create basic directories, copy raw read files, perform initial FastQC analysis, 
# trimming and final FastQC analysis
./trimming.sh
#================================================================================#


#================================================================================#
# run READemption pipeline, this includes: mapping to reference genome, gene quantification,
# generation of coverage files (full coverage and Xprime coverage), basic visualisation
# and basic DESeq2 analyis
# merged read files are generated and mapping as well as coverage files are computed 
# (full coverage and Xprime coverage)
./reademption_complete.sh
#================================================================================#


#================================================================================#
# convert raw coverage files to bedfiles using bedops
#export dir="READemption_PNPase"
./wig2bed.sh
#================================================================================#


#================================================================================#
# perform Xprime end quantification in R
echo "Xprime end quantification started..."
Rscript Xprime_quantification.R $strain1_rep1 $strain1_rep2 $strain1_rep3 $strain2_rep1 $strain2_rep2 $strain2_rep3 $cutoff_counts $path_genomic_fasta
Rscript Xprime_quantification_full.R $strain1_rep1 $strain1_rep2 $strain1_rep3 $strain2_rep1 $strain2_rep2 $strain2_rep3 $cutoff_counts $path_genomic_fasta $cutoff_ratio
echo "quantification finished"
#================================================================================#


#================================================================================#
# perform DESeq2 analysis to compute nucleotide wise fold changes
echo "Xprime DESeq analysis started..."
Rscript Xprime_DESeq.R $strain1_rep1 $strain1_rep2 $strain1_rep3 $strain2_rep1 $strain2_rep2 $strain2_rep3 $strain1 $strain2 $cutoff_log2FC $cutoff_padj $dir
#================================================================================#


#================================================================================#
# BEDtools operations:
# sort filtered output from previous DESeq2 analysis
# merge all X' ends  within a distance of 3 nucleotides
# compute intersection with annotation file
./bedtools.sh
#================================================================================#

#================================================================================#
# statistic in R
echo "R statistic analayis started..."
Rscript statistic.R $strain1 $strain2
#================================================================================#

