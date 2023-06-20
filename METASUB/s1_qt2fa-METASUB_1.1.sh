#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -o $HOME/log/
#$ -e $HOME/log/
#$ -pe def_slot 8

#==============================================================================
# qt2fa-AMR, METASUB version / created by Ryo Honda (END-AMR-Asia), 2023-06-07 
# The shell script for NIG supercomputer (Grid Engine) 
# with Singularity package of BioContainer
#==============================================================================
# This shell script creates fasta files for ARG database search from shotgun sequences by: 
# 1. quality trimming by fastp, 
# 2. merging paired-end sequences into one fasta file (without matching pairs).
#------------------------------------------------------------------------------

# ** specify the directory of other users if you need to refer.
# **他のユーザのファイルを参照する場合は次で指定 **必要ない場合は削除**
export SINGULARITY_BINDPATH="/home/user/METASUB"

####### Parameter setting ######################################################
# Your working directory 自分の作業ディレクトリ
DIR_WORKING="/home/user/METASUB"

## Listing sequence data files 
## Case 1: to input the list from a text file (1 file in 1 line).
## 配列リストをファイルから読み込む（1行1ファイルで記述）
## 【重要】改行コードはLF(UNIX/Mac)で。Windows標準ではCR+LFになってしまうので，エディタでLF(UNIX/Mac)を必ず指定。
#FILE_LIST="${DIR_WORKING}/seqlist.txt"
#IFS=$'\n'  # CR+LF(Windows)の場合は '\r\n' にする
#LIST=(`cat ${FILE_LIST}`)

## Case 2: for sequential numbered files 連番のついたファイルの場合
## Sample name prefix of sequence files, 配列ファイル名の最初についているサンプル名
PREFIX="DNBSEQ_"
## Numbering of samples to be processed 処理する配列ファイルの番号
START=100
END=101
for ((i=0;i<=`expr $END - $START`;i++)); do
	n=`expr $i + $START`; LIST[$i]="${PREFIX}${n}"
done

# The directory of raw sequence files (specify the absolute path. do not include the final '/') 
# 配列生データのあるディレクトリ（絶対パスで指定。最後にスラッシュ '/' は含めない）
DIR_RAW="/home/user/METASUB/0.raw"

# The directories for output files (specify the absolute path. do not include the final '/')  
# 配列ファイルを出力するディレクトリ（絶対パスで指定。最後にスラッシュ '/' は含めない）
DIR_QT="${DIR_WORKING}/1.trimmed"
DIR_FA="${DIR_WORKING}/2.fasta"

# location of the container image for singularity 
FASTP="/usr/local/biotools/f/fastp:0.23.3--h5f740d0_0"

####### Run program #############################################################
# Execute from Singularity. Singularityから実行
# Check the file names in the command line below are correct. 

DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
i=0; n=${#LIST[@]}
echo "[${DATE}] $JOB_NAME started. (${i}/${n})"
# make output directories if not existed
mkdir -p ${DIR_QT} ${DIR_FA}
dir_tmp="${DIR_FA}/tmp"; mkdir -p ${dir_tmp}; trap 'rm -R ${dir_tmp}' 0

for SAMPLE in "${LIST[@]}"; do
	# quality trimming by fastp
	singularity exec ${FASTP} fastp\
	--in1 ${DIR_RAW}/${SAMPLE}_Read1.fq.gz\
	--in2 ${DIR_RAW}/${SAMPLE}_Read2.fq.gz\
	--out1 ${DIR_QT}/${SAMPLE}_R1.qt.fq\
	--out2 ${DIR_QT}/${SAMPLE}_R2.qt.fq\
	--html ${DIR_QT}/${SAMPLE}.report.html\
	--json ${DIR_QT}/${SAMPLE}.report.json\
	-q 20 -t 1 -T 1 -l 20 -w 8
	
	# convert fastq to fasta
	awk '(NR - 1) % 4 < 2' ${DIR_QT}/${SAMPLE}_R1.qt.fq | sed 's/@/>/' > ${dir_tmp}/${SAMPLE}_R1.fa
	awk '(NR - 1) % 4 < 2' ${DIR_QT}/${SAMPLE}_R2.qt.fq | sed 's/@/>/' > ${dir_tmp}/${SAMPLE}_R2.fa
	# interpose two fasta files (R1&R2) into one fasta file (R12).
	cat ${dir_tmp}/${SAMPLE}_R1.fa ${dir_tmp}/${SAMPLE}_R2.fa > ${DIR_FA}/${SAMPLE}_R12.fa
	rm ${dir_tmp}/${SAMPLE}_R1.fa ${dir_tmp}/${SAMPLE}_R2.fa
 	# gzip the filtered fastq file
  	gzip ${DIR_QT}/${SAMPLE}_R1.qt.fq ${DIR_QT}/${SAMPLE}_R2.qt.fq
	
	# record the progress
	DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
	i=`expr $i + 1`
	echo "[${DATE}] ${SAMPLE} is done. (${i}/${n})"
done
# record the completion of the script
DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
echo "[$DATE] $JOB_NAME completed. (${i}/${n})"
exit 0

