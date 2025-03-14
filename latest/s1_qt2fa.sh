#!/bin/bash
JOB_NAME=$SLURM_JOB_NAME # only for Slurm job scheduler
#==============================================================================
# qt2fa-AMR ver.2 / created by Ryo Honda, Last updated: 2025-03-14
# The shell script for NIG supercomputer with the Singularity package of BioContainer
#==============================================================================
# This shell script creates fasta files for ARG database search from shotgun sequences by: 
# 1. quality trimming by fastp, 
# 2. merging paired-end sequences into one fasta file (without matching pairs).
#------------------------------------------------------------------------------

threads=8 # CPU threads / スレッド数
####### Parameter setting ######################################################
# Your working directory 作業ディレクトリ（配列データが格納されている場所）
DIR_WORKING="/home/ryohonda/GlobalAMR"

# Directories and files
## Suffix of raw sequence files / Raw配列ファイルの語尾
sfx_r1="_1.fq" # suffix of read 1 of the paired-end sequence
sfx_r2="_2.fq" # suffix of read 2 of the paired-end sequence

## The directory of raw sequence files (specify the absolute path. do not include the final '/') 
## 配列生データのあるディレクトリ（絶対パスで指定。最後にスラッシュ '/' は含めない）
DIR_RAW="${DIR_WORKING}/0.raw"
## The directory to output filtered sequence files
## トリミング後の配列データを出力するディレクトリ（絶対パスで指定。最後にスラッシュ '/' は含めない）
DIR_QT="${DIR_WORKING}/1.qt"
## The output directory of fasta files for blast (specify the absolute path. do not include the final '/')  
## Blast用のFastaファイルを出力するディレクトリ（絶対パスで指定。最後にスラッシュ '/' は含めない）
DIR_FA="${DIR_WORKING}/2.fasta"

# Sample list
## Choose listing method of sequence data files (choose 0 or 1)
## 配列ファイルリストの指定方法 (0か1を選択）
LISTING="0" 
#----------------------------------------
# 0: to use all the raw files in the directory specified as DIR_RAW above.
#	 DIR_RAWで指定したディレクトリーのすべてのファイルを使用する。
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

#====== Singularity (AppContainer) settings ==============================
## Singularity: ** specify the directory of other users if you need to refer.
## Singularity: **他のユーザのファイルを参照する場合は次で指定 **必要ない場合は削除**
export SINGULARITY_BINDPATH="${DIR_WORKING}"
## Singularity: location of the container image of the singularity package
FASTP="/usr/local/biotools/f/fastp:0.23.4--h125f33a_4"

####### Run program #############################################################
## Listing the raw sequence files.
if [ $LISTING -eq 0 ]; then 
	# Case 0: get the sample list from the fastq file list of the directory 
	LIST=($(ls ${DIR_RAW}/*${sfx_r1} ${DIR_RAW}/*${sfx_r2} | sed -E "s/(${sfx_r1}|${sfx_r2})$//" | \
	xargs -n 1 basename | sort -u))
else
	# Case 1: get the list from a file. 
	LIST=(`cat ${FILE_LIST}`)
fi
## Execute from Singularity. Singularityから実行
## Check the file names in the command line below are correct. 
DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
i=0; n=${#LIST[@]}
echo "[${DATE}] $JOB_NAME started. (${i}/${n})"

## make output directories if not existed
mkdir -p ${DIR_QT} ${DIR_FA}
dir_tmp="${DIR_FA}/tmp"; mkdir -p ${dir_tmp}; trap 'rm -R ${dir_tmp}' 0

for SAMPLE in "${LIST[@]}"; do
	## quality trimming by fastp
	singularity exec ${FASTP} \
	fastp\
	--in1 ${DIR_RAW}/${SAMPLE}${sfx_r1}\
	--in2 ${DIR_RAW}/${SAMPLE}${sfx_r2}\
	--out1 ${DIR_QT}/${SAMPLE}_R1.qt.fq\
	--out2 ${DIR_QT}/${SAMPLE}_R2.qt.fq\
	--html ${DIR_QT}/${SAMPLE}.report.html\
	--json ${DIR_QT}/${SAMPLE}.report.json\
	-q 20 -t 1 -T 1 -l 20 -w ${threads}
	## convert fastq to fasta
	awk '(NR - 1) % 4 < 2' ${DIR_QT}/${SAMPLE}_R1.qt.fq | sed 's/@/>/' > ${dir_tmp}/${SAMPLE}_R1.fa
	awk '(NR - 1) % 4 < 2' ${DIR_QT}/${SAMPLE}_R2.qt.fq | sed 's/@/>/' > ${dir_tmp}/${SAMPLE}_R2.fa
	## concatenate two fasta files (R1&R2) into one fasta file (R12).
	cat ${dir_tmp}/${SAMPLE}_R1.fa ${dir_tmp}/${SAMPLE}_R2.fa > ${DIR_FA}/${SAMPLE}_R12.fa
	#gzip the filtered fastq file
	#singularity exec /usr/local/biotools/p/pigz:2.3.4\
	#pigz -p ${threads} ${DIR_QT}/${SAMPLE}_R1.qt.fq ${DIR_QT}/${SAMPLE}_R2.qt.fq 
	
	## record the progress
	DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
	i=`expr $i + 1`
	echo "[${DATE}] ${SAMPLE} is done. (${i}/${n})"
done
## record the completion of the script
DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
echo "[$DATE] $JOB_NAME completed. (${i}/${n})"
exit 0

 
