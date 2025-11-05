#!/bin/bash

# To calculate well counts for Fig2B

## To run script
	## parallel -j8 "bash fastq_to_wellcounts.sh" {} :::: references/igvf-sequence-files.txt
	## Input references/igvf-sequence-files.txt is a file of all the fastq IGVF filepaths.  

# Export function to be used by GNU parallel
gs_path="$1"
gs_path=$(echo "$gs_path" | tr -d '\r' | xargs)
filename=$(basename "$gs_path")
outfile="../well-counts/${filename%.fastq.gz}.counts.txt"

echo "Processing $filename..."

curl -L -o - "$gs_path" | \
gzcat - | \
awk 'NR % 4 == 2' | \
cut -c116-123 | \
rev | tr 'ACGT' 'TGCA' | \
sort | uniq -c > "$outfile"

echo "Done $outfile"