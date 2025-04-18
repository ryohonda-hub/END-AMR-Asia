#!/bin/bash
JOB_NAME="${SLURM_JOB_NAME}" # any job name for your reference
#========================================
# Convert fastq files to concatenated fasta files / created by Ryo Honda, Last updated: 2025-03-14
#
# 1. All fastq files in a specified directory will be converted to fasta format
# 2. Paired-end fasta files are concatenated into one fasta file.
#========================================

# Set the working directory                                                                                                                           
DIR_WORKING="/home/ryohonda/GlobalAMR"

## Suffix of raw sequence files / Raw配列ファイルの語尾
sfx_r1="_R1.qt.fq" # suffix of read 1 of the paired-end sequence
sfx_r2="_R2.qt.fq" # suffix of read 2 of the paired-end sequence

# Input and output directories
DIR_QT="${DIR_WORKING}/1.qt"
DIR_FA="${DIR_WORKING}/2.fasta"

# Sample list
## Choose listing method of sequence data files (choose 0 or 1)
## 配列ファイルリストの指定方法 (0か1を選択）
LISTING="0" 
#----------------------------------------
# 0: to use all the quality-trimmed files in the directory specified as DIR_RAW above.
#	 DIR_QTで指定したディレクトリーのすべてのファイルを使用する。
#----------------------------------------
# 1: to input the file list from a text file (1 file in 1 line).
#    配列リストをファイルから読み込む（1行1ファイルで記述）
#
## For Case 1: Specify the text file of sequence list and newline code
## 1の場合: 配列名一覧のファイルと改行コードを指定
FILE_LIST="${DIR_WORKING}/sralist-all.txt"
IFS=$'\n' # '\n' for Mac/Unix, '\r\n' for Windows.
## [IMPORTANT] Specify the correct newline code in IFS. 
## 【重要】ファイルに使われている正しい改行コードをIFSに指定。
##  UNIX/Mac uses LF(\n); Windows uses CR+LF (\r\n)
#----------------------------------------

####### Run program #############################################################
## Listing the raw sequence files.
if [ $LISTING -eq 0 ]; then 
	# Case 0: get the sample list from the fastq file list of the directory 
	LIST=($(ls ${DIR_QT}/*${sfx_r1} ${DIR_QT}/*${sfx_r2} | sed -E "s/(${sfx_r1}|${sfx_r2})$//" | \
	xargs -n 1 basename | sort -u))
else
	# Case 1: get the list from a file. 
	LIST=(`cat ${FILE_LIST}`)
fi
echo "SAMPLE LIST:"; echo "${LIST[@]}"

## record the start of the script
DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
i=0; n=${#LIST[@]}
echo "[${DATE}] $JOB_NAME started. (${i}/${n})"

## make output directories if not existed
mkdir -p ${DIR_FA}
dir_tmp="${DIR_FA}/tmp"; mkdir -p ${dir_tmp}; trap 'rm -R ${dir_tmp}' 0

# Convert all fastq files into fasta format and save them in the tmp folder.
for SAMPLE in "${LIST[@]}"; do
	## convert fastq to fasta
	awk '(NR - 1) % 4 < 2' ${DIR_QT}/${SAMPLE}_R1.qt.fq | sed 's/@/>/' > ${dir_tmp}/${SAMPLE}_R1.fa
	awk '(NR - 1) % 4 < 2' ${DIR_QT}/${SAMPLE}_R2.qt.fq | sed 's/@/>/' > ${dir_tmp}/${SAMPLE}_R2.fa
	## concatenate two fasta files (R1&R2) into one fasta file (R12).
	cat ${dir_tmp}/${SAMPLE}_R1.fa ${dir_tmp}/${SAMPLE}_R2.fa > ${DIR_FA}/${SAMPLE}_R12.fa
	## record the progress
	DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
	i=`expr $i + 1`
	echo "[${DATE}] ${SAMPLE} is done. (${i}/${n})"
done
## record the completion of the script
DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
echo "[$DATE] $JOB_NAME completed. (${i}/${n})"
exit 0
