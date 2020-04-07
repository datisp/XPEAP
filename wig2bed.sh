#!/bin/bash

#conversion wig2bed

echo "script_wig2bed.sh started..."
date
echo "..."

dir_out=3prime_analysis/wig2bed

source activate $conda_bedops

cd $dir/output/coverage_3prime/coverage-raw
for f in *.wig; do
	wig2bed < $f > ../../../../$dir_out/$f".bed";
done

echo "script_wig2bed.sh finished."
date
