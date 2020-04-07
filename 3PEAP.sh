#!/bin/bash

echo " ________                        .__                                              .__               .__        "
echo " \_____  \          _____________|__| _____   ____           _____    ____ _____  |  | ___.__. _____|__| ______"
echo "   _(__  <   ______ \____ \_  __ \  |/     \_/ __ \   ______ \__  \  /    \\__  \ |  |<   |  |/  ___/  |/  ___/"
echo "  /       \ /_____/ |  |_> >  | \/  |  Y Y  \  ___/  /_____/  / __ \|   |  \/ __ \|  |_\___  |\___ \|  |\___ \ "
echo " /______  /         |   __/|__|  |__|__|_|  /\___  >         (____  /___|  (____  /____/ ____/____  >__/____  >"
echo "        \/          |__|                  \/     \/               \/     \/     \/     \/         \/        \/"

# 2020-03-17 version 0.1
# written by Daniel-Timon Spanka
# contact daniel-timon.spanka@mikro.bio.uni-giessen.de


#================================================================================#
# requirements:
# FastQC >= v0.11.8
# MultiQC >= v0.9 in conda environment
# READemption >= v0.4.3 in conda environment
# Bedops >= v2.4.37 in conda environment
# BEDtools >= v2.25.0
#
# requirements for R:
# tidyverse
# DESeq2
# gplots
# RColorBrewer
#================================================================================#


#================================================================================#
# define all global variables
# path to raw data in fastq file format
export path_raw_data="/nfs/agmikmol/Timon/raw_reads/2019-07-02_RNAseq_raw_wt_pnp/2018-09-27_FastQ_P180084"

# path to gff and genomic.fasta
export path_gff="/nfs/agmikmol/Timon/refseq/GCF_000012905.2_ASM1290v2_genomic_extended.gff"
export path_genomic_fasta="/nfs/agmikmol/Timon/refseq/GCF_000012905.2_ASM1290v2_genomic.fa"

# path to fastqc bin file
export path_fastqc="/homes/dtspanka/Tools/FastQC"

# unique regular expression pattern for replicates of both conditions, adjust according to filenames
export var_strain1="pnp"
export var_strain2="WT_neu"

# variables for R analysis
# save filenames of each read file (without file type extension!)
export strain1_rep1="L1802003_Nr25_pnp_exp_1"
export strain1_rep2="L1802004_Nr26_pnp_exp_2"
export strain1_rep3="L1802005_Nr27_pnp_exp_3"
export strain2_rep1="L1802021_Nr43_WT_neu_1"
export strain2_rep2="L1802022_Nr44_WT_neu_2"
export strain2_rep3="L1802023_Nr45_WT_neu_3"
export cutoff_counts="10"
export cutoff_padj="0.05"
export cutoff_log2FC="1"
export condition1="pnp"
export condition2="wt"

# number of cores to be used for computation
export n_cores="12"

# names of conda environments
export conda_multiqc="multiqc"
export conda_trim_galore="trimming"
export conda_reademption="rnaseq"
export conda_bedops="bedops"

# name for READemption project folders
export dir="READemption_PNPase"

# READemption settings
#export input_libs_old="L1802021_Nr43_WT_neu_1_trimmed.fq,L1802022_Nr44_WT_neu_2_trimmed.fq,L1802023_Nr45_WT_neu_3_trimmed.fq,L1802003_Nr25_pnp_exp_1_trimmed.fq,L1802004_Nr26_pnp_exp_2_trimmed.fq,L1802005_Nr27_pnp_exp_3_trimmed.fq"
export input_libs=$strain2_rep1"_trimmed.fq,"$strain2_rep2"_trimmed.fq,"$strain2_rep3"_trimmed.fq,"$strain1_rep1"_trimmed.fq,"$strain1_rep2"_trimmed.fq,"$strain1_rep3"_trimmed.fq"
export input_conditions=$condition2","$condition2","$condition2","$condition1","$condition1","$condition1
#================================================================================#


#================================================================================#
# create basic directories, copy raw read files, perform initial FastQC analysis, 
# trimming and final FastQC analysis
#./trimming.sh
#================================================================================#


#================================================================================#
# run READemption pipeline, this includes: mapping to reference genome, gene quantification,
# generation of coverage files (full coverage and 3-prime coverage), basic visualisation
# and basic DESeq2 analyis
# Merged read files are generated and mapping as well as coverage files are computeted 
# (full coverage and 3-prime coverage)
#./reademption_complete.sh
#================================================================================#


#================================================================================#
# convert raw coverage files to bedfiles using bedops
#export dir="READemption_PNPase"
#./wig2bed.sh
#================================================================================#


#================================================================================#
# perform 3prime end quantification in R
echo "3prime end quantification started..."
#Rscript 3prime_quantification.R $strain1_rep1 $strain1_rep2 $strain1_rep3 $strain2_rep1 $strain2_rep2 $strain2_rep3 $cutoff_counts $path_genomic_fasta
echo "quantification finished"
#================================================================================#


#================================================================================#
# perform DESeq2 analysis to compute nucleotide wise fold changes
echo "3prime DESeq analysis started..."
#Rscript 3prime_DESeq.R $strain1_rep1 $strain1_rep2 $strain1_rep3 $strain2_rep1 $strain2_rep2 $strain2_rep3 $condition1 $condition2 $cutoff_log2FC $cutoff_padj $dir
#================================================================================#


#================================================================================#
# BEDtools operations:
# sort filtered output from previous DESeq2 analysis
# merge all 3'end  within a distance of 3 nucleotides
# compute intersection with annotation file
#./bedtools.sh
#================================================================================#

#================================================================================#
# statistic in R
echo "R statistic analayis started..."
#Rscript statistic.R $condition1 $condition2
#================================================================================#

