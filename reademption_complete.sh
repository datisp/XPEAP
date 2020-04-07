#!/bin/bash

### READemption pipeline for PNPase versus WT new
### 2020-02-19

echo "script_reademption_complete.sh started..."

# activate READemption conda env
source activate $conda_reademption
#bash
#source activate rnaseq

reademption create $dir

date > ./$dir/log_align

# READemption only acceptes fastq -> unzip fq.gz files
# only once necessary
gunzip ./reads_trimmed/*.fq.gz


# copy .fa, gff and reads to respective input folders
#ln -s /nfs/agmikmol/Timon/PNPase/reads_trimmed/*.fq ./$dir/input/reads/.
#ln -s /nfs/agmikmol/Timon/refseq/GCF_000012905.2_ASM1290v2_genomic.fa ./$dir/input/reference_sequences/.
#ln -s /nfs/agmikmol/Timon/refseq/GCF_000012905.2_ASM1290v2_genomic_extended.gff ./$dir/input/annotations/.
cp ./reads_trimmed/*trimmed.fq ./$dir/input/reads/
cp $path_genomic_fasta ./$dir/input/reference_sequences/ref_seq.fa
cp $path_gff ./$dir/input/annotations/annotation.gff

# READemption align
reademption align -p $n_cores --progress --fastq ./$dir &>> ./$dir/log

# READemption gene quanti
date > ./$dir/log_gene_quanti
reademption gene_quanti -p $n_cores --features CDS,tRNA,rRNA,tmRNA,ncRNA,sRNA ./$dir &>> ./$dir/log_gene_quanti

# READemption
date > ./$dir/log_deseq
#input_libs=L1802003_Nr25_pnp_exp_1_trimmed.fq,L1802004_Nr26_pnp_exp_2_trimmed.fq,L1802005_Nr27_pnp_exp_3_trimmed.fq,L1802021_Nr43_WT_neu_1_trimmed.fq,L1802022_Nr44_WT_neu_2_trimmed.fq,L1802023_Nr45_WT_neu_3_trimmed.fq
#input_conditions=PNP_exp,PNP_exp,PNP_exp,WT_new_exp,WT_new_exp,WT_new_exp
reademption deseq --libs $input_libs --conditions $input_conditions ./$dir &>> ./$dir/log_deseq

# READemption visualisation
#reademption viz_align ./$dir
#reademption viz_gene_quanti ./$dir
reademption viz_deseq ./$dir
reademption coverage -p $n_cores ./$dir

# rename coverage folder
mv ./$dir/output/coverage ./$dir/output/coverage_full

# generate new folder structure for 3prime mapping
mkdir -p ./$dir/output/coverage/{coverage-raw,coverage-tnoar_mil_normalized,coverage-tnoar_min_normalized}

# generate 3prime mapping
reademption coverage -p $n_cores --coverage_style last_base_only ./$dir
mv ./$dir/output/coverage ./$dir/output/coverage_3prime


#--------------------------------------------------------------------------------#
##### 

# folder structure merged files
dir=$dir"_merged"
reademption create $dir
cp ./reads_trimmed/*trimmed.fq ./$dir/input/reads/

# merge read files
dir_reads=$dir"/input/reads"
cat ./$dir_reads/*$var_strain1*.fq > ./$dir_reads/$var_strain1"_merged.fq"
cat ./$dir_reads/*$var_strain2*.fq > ./$dir_reads/$var_strain2"_merged.fq"
rm -r ./$dir_reads/*trimmed.fq
cp $path_genomic_fasta/GCF_000012905.2_ASM1290v2_genomic.fa ./$dir/input/reference_sequences/
cp $path_gff/GCF_000012905.2_ASM1290v2_genomic_extended.gff ./$dir/input/annotations/

# READemption align
reademption align -p $n_cores --progress --fastq ./$dir &>> ./$dir/log

#READemption converage
reademption coverage -p $n_cores ./$dir

# rename coverage folder
mv ./$dir/output/coverage ./$dir/output/coverage_full_merged

# generate new folder structure for 3prime mapping
mkdir -p ./$dir/output/coverage/{coverage-raw,coverage-tnoar_mil_normalized,coverage-tnoar_min_normalized}

# generate 3prime mapping
reademption coverage -p $n_cores --coverage_style last_base_only ./$dir
mv ./$dir/output/coverage ./$dir/output/coverage_3prime_merged

echo "script_reademption_complete.sh finished."
date
