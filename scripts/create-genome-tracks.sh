#!/bin/bash

#create genome tracks

## required dependencies:
	## deeptools: pip install deeptools
	## sinto: pip install sinto
	## SparK.py: weten https://raw.githubusercontent.com/harbourlab/SparK/refs/heads/master/SparK.py
	## GNU parallel: https://ftpmirror.gnu.org/parallel/parallel-latest.tar.bz2
	## samtools: https://www.htslib.org/download


echo "Merging BAM files."
samtools merge -@16 share-v2-atac-merged.bam Subpool_1.bam Subpool_2.bam Subpool_3.bam Subpool_4.bam Subpool_5.bam 
samtools sort -@6 -m 2G -o share-v2-atac-merged-sorted.bam share-v2-atac-merged.bam 
rm share-v2-atac-merged.bam
samtools index -@8 share-v2-atac-merged-sorted.bam

echo "Finished merging BAM files."

echo "Filter merged BAM file." 
sinto filterbarcodes -b share-v2-atac-merged-sorted.bam -c results/barcode_celltype.txt -p 8 --outdir celltype_bams/
Cd celltype_bams
echo "Finished filtering merged BAM file." 

echo "Sort cell type-specific BAM." 
input_dir="."
output_dir="./sorted_bams"
mkdir -p "$output_dir"

# Loop through all BAM files in the input directory
for bam_file in "$input_dir"/*.bam; do
	base_name=$(basename "$bam_file" .bam)
	output_file="$output_dir/${base_name}.sorted.bam"
	echo "Sorting $bam_file to $output_file"
        samtools sort -m 2G -@8 -o "$output_file" "$bam_file"
	samtools index -@8 "$output_file"
done
echo "Finished sorting BAM files."

echo "Convert bam into normalized bedgraph"
ls *.bam | parallel 'bamCoverage -b {} -o {.}.norm.bdg -bs 1 -of bedgraph --normalizeUsing CPM'
ls *.bdg | parallel 'bgzip {}; tabix -p bed {}.gz; echo {}'

echo "Create genome track SVGs"

python3 SparK.py -pr chr2:102500066-102610067 -cf GABAergic_interneuron.sorted.norm.bdg.gz barrier_cell.sorted.norm.bdg.gz glutamatergic_neuron.sorted.norm.bdg.gz microglial_cell.sorted.norm.bdg.gz astrocyte.sorted.norm.bdg.gz epithelial_cell.sorted.norm.bdg.gz medium_spiny_neuron.sorted.norm.bdg.gz oligodendrocyte.sorted.norm.bdg.gz -gtf genes.gtf -gl GABAergic_interneuron barrier_cell glutamatergic_neuron microglial_cell astrocyte epithelial_cell medium_spiny_neuron oligodendrocyte -dg Slc1a2 -sm 5 -bed S1c1a2.bed S1c1a2-tss.bed -w 300 -cs 2 2 2 2 2 2 2 2 -o s1c1a2-tracks-norm

python3 SparK.py -pr chr17:10425156-10535157 -cf GABAergic_interneuron.sorted.norm.bdg.gz barrier_cell.sorted.norm.bdg.gz glutamatergic_neuron.sorted.norm.bdg.gz microglial_cell.sorted.norm.bdg.gz astrocyte.sorted.norm.bdg.gz epithelial_cell.sorted.norm.bdg.gz medium_spiny_neuron.sorted.norm.bdg.gz oligodendrocyte.sorted.norm.bdg.gz -gtf genes.gtf -gl GABAergic_interneuron barrier_cell glutamatergic_neuron microglial_cell astrocyte epithelial_cell medium_spiny_neuron oligodendrocyte -dg Qki -sm 5 -cs 1 1 1 1 1 1 1 1 -bed qki.bed qki-tss.bed -w 300 -o qki-tracks-norm
