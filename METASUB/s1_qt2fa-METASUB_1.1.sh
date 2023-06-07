#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -o $HOME/log/
#$ -e $HOME/log/
#$ -pe def_slot 8

#==============================================================================
# qt2fa-AMR, METASUB version / created by Ryo Honda (END-AMR-Asia), 2023-05-27 
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
# Sample name prefix of sequence files, 配列ファイル名の最初についているサンプル名
SAMPLE="DNBSEQ_"
# Numbering of samples to be processed 処理する配列ファイルの番号
START=270
END=270

# The directory of raw sequence files (specify the absolute path. do not include the final '/') 
# 配列生データのあるディレクトリ（絶対パスで指定。最後にスラッシュ '/' は含めない）
DIR_RAW="/home/user/METASUB/0.raw"

# The directories for output files (specify the absolute path. do not include the final '/')  
# 配列ファイルを出力するディレクトリ（絶対パスで指定。最後にスラッシュ '/' は含めない）
DIR_OUT="/home/user/METASUB"
DIR_QT="${DIR_OUT}/1.trimmed"
DIR_FA="${DIR_OUT}/2.fasta"

# location of the container image for singularity 
FASTP="/usr/local/biotools/f/fastp:0.23.3--h5f740d0_0"

####### Run program #############################################################
# Execute from Singularity. Singularityから実行
# Check the file names in the command line below are correct. 

DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
echo "[${DATE}] $JOB_NAME started. ${SAMPLE}${START}-${END}"
# make output directories if not existed
mkdir -p ${DIR_QT} ${DIR_FA}
dir_tmp="${DIR_FA}/tmp"; mkdir -p ${dir_tmp}; trap 'rm -R ${dir_tmp}' 0

for ((i=START; i<=END; i++))
do
	# quality trimming by fastp
	singularity exec ${FASTP} fastp\
	--in1 ${DIR_RAW}/${SAMPLE}${i}_Read1.fq.gz\
	--in2 ${DIR_RAW}/${SAMPLE}${i}_Read2.fq.gz\
	--out1 ${DIR_QT}/${SAMPLE}${i}_Read1.qt.fq.gz\
	--out2 ${DIR_QT}/${SAMPLE}${i}_Read2.qt.fq.gz\
	--html ${DIR_QT}/${SAMPLE}${i}.report.html\
	--json ${DIR_QT}/${SAMPLE}${i}.report.json\
	-q 20 -t 1 -T 1 -l 20 -w 8
	
	# convert fastq to fasta
	gunzip -c ${DIR_QT}/${SAMPLE}${i}_Read1.qt.fq.gz | awk '(NR - 1) % 4 < 2' | sed 's/@/>/' > ${dir_tmp}/${SAMPLE}${i}_R1.fa
	gunzip -c ${DIR_QT}/${SAMPLE}${i}_Read2.qt.fq.gz | awk '(NR - 1) % 4 < 2' | sed 's/@/>/' > ${dir_tmp}/${SAMPLE}${i}_R2.fa
	# interpose two fasta files (R1&R2) into one fasta file (R12).
	cat ${dir_tmp}/${SAMPLE}${i}_R1.fa ${dir_tmp}/${SAMPLE}${i}_R2.fa > ${DIR_FA}/${SAMPLE}${i}_R12.fa
	rm ${dir_tmp}/${SAMPLE}${i}_R1.fa ${dir_tmp}/${SAMPLE}${i}_R2.fa
	
	# record the progress
	DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
	echo "[${DATE}] ${SAMPLE}${i} is done."
done
# record the completion of the script
DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
echo "[$DATE] $JOB_NAME completed."
exit 0

