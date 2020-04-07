#!/bin/bash

### processing, trimming, quality control
### 2020-02-19

echo "script_trimming.sh started..."
date

# create folder structure 
mkdir reads_raw
mkdir reads_trimmed
mkdir -p 3prime_analysis/{3prime_DESeq,wig2bed}

# copy raw reads
echo "copying raw data..."
cp $path_raw_data/*$var_strain1*fq.gz ./reads_raw/
cp $path_raw_data/*$var_strain2*fq.gz ./reads_raw/

# perform FastQC analysis for unprocessed reads
echo "performing FastQC analysis for raw reads..."
$path_fastqc/fastqc ./reads_raw/*.fq.gz

# activate conda for MultiQC analysis
source activate $conda_multiqc

# perform MultiQC analysis for raw reads
multiqc ./reads_raw/*fastqc.zip -o ./reads_raw/

# activate conda trimming for trim-galore
echo "removal of adapter sequences and quality trimming..." 
source activate $conda_trim_galore
trim_galore -o ./reads_trimmed -q 10 --length 15 -j 8 ./reads_raw/*.gz &> ./reads_trimmed/log

# perform FastQC for trimmed reads
echo "performing FastQC analysis for trimmed reads..."
$path_fastqc/fastqc ./reads_trimmed/*fq.gz

# activate conda for MultiQC analysis
source activate $conda_multiqc

# perform MultiQC analysis for raw reads
multiqc ./reads_trimmed/*fastqc.zip -o ./reads_trimmed/

# READemption only allows unzipped fastq-files
gunzip ./reads_trimmed/*.fq.gz

echo "script_trimming.sh finished."
date
