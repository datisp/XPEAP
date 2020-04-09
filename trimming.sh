#!/bin/bash

### processing, trimming, quality control
### 2020-04-07

echo "script_trimming.sh started..."
date

# create folder structure 
mkdir reads_raw
mkdir reads_trimmed

# copy raw reads
echo "copying raw data..."
#cp $path_raw_data/*$var_strain1*fq.gz ./reads_raw/
#cp $path_raw_data/*$var_strain2*fq.gz ./reads_raw/

cp $path_raw_data/{$strain1_rep1$filetype$compression,$strain1_rep2$filetype$compression,$strain1_rep3$filetype$compression,$strain2_rep1$filetype$compression,$strain2_rep2$filetype$compression,$strain2_rep3$filetype$compression} ./reads_raw/

# perform FastQC analysis for unprocessed reads
echo "performing FastQC analysis for raw reads..."
$path_fastqc/fastqc ./reads_raw/*$filetype$compression

# activate conda for MultiQC analysis
source activate $conda_multiqc

# perform MultiQC analysis for raw reads
multiqc ./reads_raw/*fastqc.zip -o ./reads_raw/

# activate conda trimming for trim-galore
echo "removal of adapter sequences and quality trimming..." 
source activate $conda_trim_galore
trim_galore -o ./reads_trimmed -q 10 --length 15 -j 8 ./reads_raw/*$compression &> ./reads_trimmed/log

# perform FastQC for trimmed reads
echo "performing FastQC analysis for trimmed reads..."
$path_fastqc/fastqc ./reads_trimmed/*$filetype$compression

# activate conda for MultiQC analysis
source activate $conda_multiqc

# perform MultiQC analysis for raw reads
multiqc ./reads_trimmed/*fastqc.zip -o ./reads_trimmed/

# READemption only allows unzipped fastq-files
gunzip ./reads_trimmed/*$filetype$compression

echo "script_trimming.sh finished."
date
