#!/bin/bash

echo "     ___          ___       ___          ___          ___   "
echo "    /__/|        /  /\     /  /\        /  /\        /  /\  "
echo "   |  |:|       /  /::\   /  /:/_      /  /::\      /  /::\ "
echo "   |  |:|      /  /:/\:\ /  /:/ /\    /  /:/\:\    /  /:/\:\"
echo " __|__|:|     /  /:/~/://  /:/ /:/_  /  /:/~/::\  /  /:/~/:/"
echo "/__/::::\____/__/:/ /://__/:/ /:/ /\/__/:/ /:/\:\/__/:/ /:/ "
echo "   ~\~~\::::/\  \:\/:/ \  \:\/:/ /:/\  \:\/:/__\/\  \:\/:/  "
echo "    |~~|:|~~  \  \::/   \  \::/ /:/  \  \::/      \  \::/   "
echo "    |  |:|     \  \:\    \  \:\/:/    \  \:\       \  \:\   "
echo "    |  |:|      \  \:\    \  \::/      \  \:\       \  \:\  "
echo "    |__|/        \__\/     \__\/        \__\/        \__\/  "

# 2020-04-09 version 1.0

#================================================================================#

# please modify the content of the following variables:

# path to folder with raw data
export path_raw_data="/nfs/agmikmol/Timon/raw_reads/2019-07-02_RNAseq_raw_wt_pnp/2018-09-27_FastQ_P180084"

# path to gff and genomic.fasta
export path_gff="/nfs/agmikmol/Timon/refseq/GCF_000012905.2_ASM1290v2_genomic_extended.gff"
export path_genomic_fasta="/nfs/agmikmol/Timon/refseq/GCF_000012905.2_ASM1290v2_genomic.fa"

export strain1="pnp"
export strain2="wt"

# save filenames of each read file (without file type extension!)
export strain1_rep1="L1802003_Nr25_pnp_exp_1"
export strain1_rep2="L1802004_Nr26_pnp_exp_2"
export strain1_rep3="L1802005_Nr27_pnp_exp_3"
export strain2_rep1="L1802021_Nr43_WT_neu_1"
export strain2_rep2="L1802022_Nr44_WT_neu_2"
export strain2_rep3="L1802023_Nr45_WT_neu_3"
export filetype=".fq"
export compression=".gz"

# path to fastqc bin file
export path_fastqc="/homes/dtspanka/Tools/FastQC"

# names of conda environments
export conda_multiqc="multiqc"
export conda_trim_galore="trimming"
export conda_reademption="rnaseq"
export conda_bedops="bedops"

export cutoff_counts="10"
export cutoff_padj="0.05"
export cutoff_log2FC="1"

# READemption settings:
export n_cores="24" # number of cores to be used for computation
export dir="READemption_PNPase" # name for READemption project folders
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
# merge all X'end  within a distance of 3 nucleotides
# compute intersection with annotation file
./bedtools.sh
#================================================================================#

#================================================================================#
# statistic in R
echo "R statistic analayis started..."
Rscript statistic.R $strain1 $strain2
#================================================================================#

