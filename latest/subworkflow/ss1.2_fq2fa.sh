#!/bin/bash
#========================================
# Convert fastq files to concatenated fasta files
#
# 1. All fastq files in a specified directory will be converted to fasta format
# 2. Paired-end fasta files are concatenated into one fasta file.
#========================================
DIR_WORKING="/home/ryohonda/JST-MIRAI"

## Suffix of raw sequence files / Raw配列ファイルの語尾
sfx_fq=".qt.fq"
sfx_r1="_R1.fa" # suffix of read 1 of the paired-end sequence
sfx_r2="_R2.fa" # suffix of read 2 of the paired-end sequence

# Input and output directories
DIR_QT="${DIR_WORKING}/1.qt"
DIR_FA="${DIR_WORKING}/2.fasta"

## make output directories if not existed
mkdir -p ${DIR_FA}
dir_tmp="${DIR_FA}/tmp"; mkdir -p ${dir_tmp}; trap 'rm -R ${dir_tmp}' 0

# Convert all fastq files into fasta format and save them in the tmp folder.
for file in ${DIR_QT}/*${sfx_fq}; do
    output_file="${dir_tmp}/${file%sfx_fq}.fa"  # 出力ファイル名の作成
    awk '(NR - 1) % 4 < 2' "$file" | sed 's/@/>/' > "$output_file"
done
# Merge paired-end files
for file in ${dir_tmp}/*.fa; do
    SAMPLE="${dir_tmp}/${file%.fa}"  # 出力ファイル名の作成
	cat ${dir_tmp}/${SAMPLE}${sfx_r1} ${dir_tmp}/${SAMPLE}${sfx_r2} > ${DIR_FA}/${SAMPLE}_R12.fa
	rm ${dir_tmp}/${SAMPLE}${sfx_r1} ${dir_tmp}/${SAMPLE}${sfx_r2}
done
exit 0
