#!/bin/bash

"""
This script takes in raw sequencing reads fastq files, grep chimeric reads that contain 
Tol2 sequences, trim off Tn5 adapter and Tol2 sequences, output the remaining flanking sequences uniquely mapped to genome (bigwig format).
Input: 
	sample fastq file
Ouput: 
	sample.bw
Example usage:
	"bash alignment_based_TransTag.sh sample.fastq.gz" 
"""

#reference genome assembly 
#zv11
REFERENCEINDEX="Genomes/zebrafish/bowtie2_index/danRer11"

#take input fastq file from command line input
ifile=$1
echo "processing file: $1"
#extract sample name	
iname=$(echo "$ifile" | awk -F'/' '{print $NF}' | sed -r 's/.fq.gz//g' | sed -r 's/.fastq.gz//g' | sed -r 's/.fastq//g' | sed -r 's/.fq//g')
echo "sample name: $iname"

#grep out chimeric reads that contain Tol2 sequences
echo "grep chimeric reads..."
less $ifile | grep -B 1 "TTTCACTTGAGTAAAATTTTTGAGTACTTTTTACACCTCTG" | grep "@" | awk '{print $1}' > "${iname}".chimeric_reads_name.txt

echo "done grep chimeric reads..."
awk 'FNR==NR {lines[$1];next} $1 in lines {c=4}c-->0' "${iname}".chimeric_reads_name.txt "${iname}"_R1.fastq > "${iname}".chimeric_reads_name.R1.fastq
awk 'FNR==NR {lines[$1];next} $1 in lines {c=4}c-->0' "${iname}".chimeric_reads_name.txt "${iname}"_R2.fastq > "${iname}".chimeric_reads_name.R2.fastq

gzip "${iname}".chimeric_reads_name.R1.fastq
gzip "${iname}".chimeric_reads_name.R2.fastq

#remove Tn5 adapter and Tol2 sequences to get flanking genomic regions
conda activate Cutadapt_env
cutadapt -b CTGTCTCT -B CTGTCTCT -o "${iname}".chimeric_reads_name.trimTn5.R1.fastq.gz -p "${iname}".chimeric_reads_name.trimTn5.R2.fastq.gz "${iname}".chimeric_reads_name.R1.fastq.gz "${iname}".chimeric_reads_name.R2.fastq.gz
cutadapt -g CTTTTTACACCTCTG -G CTTTTTACACCTCTG -o "${iname}".chimeric_reads_name.trimTol2.R1.fastq.gz -p "${iname}".chimeric_reads_name.trimTol2.R2.fastq.gz "${iname}".chimeric_reads_name.trimTn5.R1.fastq.gz "${iname}".chimeric_reads_name.trimTn5.R2.fastq.gz
conda deactivate


#map the reads and only use uniquely mapped reads to generate bigwig files
conda activate Bowtie2_mapping_env
echo "working on mapping of sample: ${iname}"
bowtie2 -p 8 -X 2000 -x $REFERENCEINDEX -1 "${iname}".chimeric_reads_name.trimTol2.R1.fastq.gz -2 "${iname}".chimeric_reads_name.trimTol2.R1.fastq.gz | samtools view -b -h -F 4 -@ 8 - | samtools sort -@ 8 - -o "${iname}".chimeric_reads_name.trimdouble.sorted.bam
samtools index "${iname}".chimeric_reads_name.trimdouble.sorted.bam

conda activate Deeptools_env
echo "working on bigwig coverage of sample: ${iname}"
bamcoverage -b "${iname}".chimeric_reads_name.trimdouble.sorted.bam -o "${iname}".chimeric_reads_name.trimdouble.q30.bw -p 8 --normalizeUsing RPKM --binSize 10 --minMappingQuality 30 
bamcoverage -b "${iname}".chimeric_reads_name.trimdouble.sorted.bam -o "${iname}".chimeric_reads_name.trimdouble.q40.bw -p 8 --normalizeUsing RPKM --binSize 10 --minMappingQuality 40
bamcoverage -b "${iname}".chimeric_reads_name.trimdouble.sorted.bam -o "${iname}".chimeric_reads_name.trimdouble.bw -p 8 --normalizeUsing RPKM --binSize 10 












