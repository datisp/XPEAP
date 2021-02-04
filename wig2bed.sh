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
cd $dir/output/coverage_Xprime/coverage-tnoar_min_normalized
for f in *.wig; do
	set -x
	mv "$f" "$(echo $f | awk -F '[_]' '{OFS = "_"}{print $1, $2, $3, $4, $5, $6, $13}')"
	set +x
done

for f in *.wig; do
	wig2bed < $f > ../../../../$dir_out/coverage_Xprime/$f".bed";
done

cd ../..

# transform full coverage
cd coverage_full/coverage-tnoar_min_normalized

for f in *.wig; do
	set -x
	mv "$f" "$(echo $f | awk -F '[_]' '{OFS = "_"}{print $1, $2, $3, $4, $5, $6, $13}')"
	set +x
done

for f in *.wig; do
	wig2bed < $f > ../../../../$dir_out/coverage_full/$f".bed";
done

echo "wig2bed.sh finished."
date
