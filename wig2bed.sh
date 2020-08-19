#!/bin/bash

#conversion wig2bed

echo "wig2bed.sh started..."
date
echo "..."

mkdir -p Xprime_analysis/{Xprime_DESeq,wig2bed}
mkdir -p Xprime_analysis/wig2bed/{coverage_Xprime,coverage_full}

dir_out=Xprime_analysis/wig2bed

source activate $conda_bedops

# transform Xprime coverage
cd $dir/output/coverage_Xprime/coverage-raw
for f in *.wig; do
	wig2bed < $f > ../../../../$dir_out/coverage_Xprime/$f".bed";
done

cd ../..

# transform full coverage
cd coverage_full/coverage-raw
for f in *.wig; do
	wig2bed < $f > ../../../../$dir_out/coverage_full/$f".bed";
done

echo "wig2bed.sh finished."
date
