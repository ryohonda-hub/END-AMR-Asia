#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -o $HOME/log/
#$ -e $HOME/log/
#$ -pe def_slot 8

#==============================================================================
# qt2fa-AMR, JST-MIRAI version / created by Ryo Honda, 2023-06-21
# The shell script for NIG supercomputer (Grid Engine) 
# with Singularity package of BioContainer
#==============================================================================
# This shell script creates fasta files for ARG database search from shotgun sequences by: 
# 1. quality trimming by fastp, 
# 2. merging paired-end sequences into one fasta file (without matching pairs).
#------------------------------------------------------------------------------

## ** specify the directory of other users if you need to refer.
## **他のユーザのファイルを参照する場合は次で指定 **必要ない場合は削除**
export SINGULARITY_BINDPATH="/home/user/JST-MIRAI"

####### Parameter setting ######################################################
# Your working directory 自分の作業ディレクトリ
DIR_WORKING="/home/user/JST-MIRAI"

#----------------------------------------
## Choose listing method of sequence data files (choose 1 or 2)
## 配列ファイルリストの指定方法 (1か2を選択）
LISTING="1" 
# 1: to input the list from a text file (1 file in 1 line).
#    配列リストをファイルから読み込む（1行1ファイルで記述）
# 2: for sequential numbered files 連番のついたファイル
#----------------------------------------
## Case 1: Specify the text file of sequence list and newline code
## 1の場合: 配列名一覧のファイルと改行コードを指定
FILE_LIST="${DIR_WORKING}/sralist.txt"
IFS=$'\n' # '\n' for Mac/Unix, '\r\n' for Windows.
## [IMPORTANT] Specify the correct newline code in IFS. 
## 【重要】ファイルに使われている正しい改行コードをIFSに指定。
#  UNIX/Mac uses LF(\n); Windows uses CR+LF (\r\n)
#----------------------------------------
## Case 2: Specify the prefix and numbering of the sequence file.
## 2の場合: 配列ファイル目の頭と，連番の開始・終了を指定
PREFIX="DRR"
START=270001
END=270005
#----------------------------------------
## Suffix of raw sequence files / Raw配列ファイルの語尾
sfx_r1="_1.fastq" # suffix of read 1 of the paired-end sequence
sfx_r2="_2.fastq" # suffix of read 2 of the paired-end sequence
#----------------------------------------

## The directory of raw sequence files (specify the absolute path. do not include the final '/') 
## 配列生データのあるディレクトリ（絶対パスで指定。最後にスラッシュ '/' は含めない）
DIR_RAW="/home/user/JST-MIRAI/0.raw"

## The directories for output files (specify the absolute path. do not include the final '/')  
## 配列ファイルを出力するディレクトリ（絶対パスで指定。最後にスラッシュ '/' は含めない）
DIR_QT="${DIR_WORKING}/1.trimmed"
DIR_FA="${DIR_WORKING}/2.fasta"

## location of the container image for singularity 
FASTP="/usr/local/biotools/f/fastp:0.23.3--h5f740d0_0"

####### Run program #############################################################
## Listing the raw sequence files.
if [ $LISTING -eq 2 ]; then 
	# in case of sequential files
	for ((i=0;i<=`expr $END - $START`;i++)); do
		n=`expr $i + $START`; LIST[$i]="${PREFIX}${n}"
	done
else
	# get the list from a file. 
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
	--in1 ${DIR_RAW}/${SAMPLE}${sfx_r1}\# 元の配列ファイルに合わせて変更
	--in2 ${DIR_RAW}/${SAMPLE}${sfx_r2}\# 元の配列ファイルに合わせて変更
	--out1 ${DIR_QT}/${SAMPLE}_R1.qt.fq\
	--out2 ${DIR_QT}/${SAMPLE}_R2.qt.fq\
	--html ${DIR_QT}/${SAMPLE}.report.html\
	--json ${DIR_QT}/${SAMPLE}.report.json\
	-q 20 -t 1 -T 1 -l 20 -w 8
	## convert fastq to fasta
	awk '(NR - 1) % 4 < 2' ${DIR_QT}/${SAMPLE}_R1.qt.fq | sed 's/@/>/' > ${dir_tmp}/${SAMPLE}_R1.fa
	awk '(NR - 1) % 4 < 2' ${DIR_QT}/${SAMPLE}_R2.qt.fq | sed 's/@/>/' > ${dir_tmp}/${SAMPLE}_R2.fa
	## interpose two fasta files (R1&R2) into one fasta file (R12).
	cat ${dir_tmp}/${SAMPLE}_R1.fa ${dir_tmp}/${SAMPLE}_R2.fa > ${DIR_FA}/${SAMPLE}_R12.fa
	#gzip the filtered fastq file
	#singularity exec /usr/local/biotools/p/pigz:2.3.4\
	#pigz -p 8 ${DIR_QT}/${SAMPLE}_R1.qt.fq ${DIR_QT}/${SAMPLE}_R2.qt.fq 
	
	## record the progress
	DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
	i=`expr $i + 1`
	echo "[${DATE}] ${SAMPLE} is done. (${i}/${n})"
done
## record the completion of the script
DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
echo "[$DATE] $JOB_NAME completed. (${i}/${n})"
exit 0

 
