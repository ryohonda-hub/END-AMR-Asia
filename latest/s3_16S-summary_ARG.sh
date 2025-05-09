#!/bin/bash
JOB_NAME="${SLURM_JOB_NAME}" # any job name for your reference
#==============================================================================
# Taxonomy analysis and sequence reads summary / created by Ryo Honda, Last updated: 2025-03-14
# The shell script for NIG supercomputer with the Singularity package of BioContainer
#==============================================================================
# This shell script creates taxonomy profile tables and a summary table of sequence reads, including:
#  - 16S-based taxonomy composition tables of each sample
#  - a table on the number of raw reads, quality reads, total 16S reads, and total ARG reads
# # This shell script recalls and requires: 
# 1. Kraken2 and Bracken2 [for taxonomy analysis using the SILVA database], 
# 2. summary_reads_16S_ARG.py [to creates a summary table of sequence reads]
#------------------------------------------------------------------------------

threads=8	# CPU threads / スレッド数
rlen=150	# sequence read length 
####### Parameter setting ######################################################
# Your working directory 自分の作業ディレクトリ
DIR_WORKING="/home/ryohonda/GlobalAMR"

# Sample list
## Choose listing method of sequence data files (choose 0 or 1)
## 配列ファイルリストの指定方法 (0か1を選択）
LISTING="0" 
#----------------------------------------
# 0: to use all the sequence files in the directory specified as DIR_QT below.
#	 DIR_QTで指定したディレクトリーのすべてのファイルを使用する。
#----------------------------------------
# 1: to input the file list from a text file (1 file in 1 line).
#    配列リストをファイルから読み込む（1行1ファイルで記述）
#
## For Case 1: Specify the text file of sequence list and newline code
## 1の場合: 配列名一覧のファイルと改行コードを指定
FILE_LIST="${DIR_WORKING}/sralist-all.txt"
## [IMPORTANT] In the sample list file,
## (1) specify 'LF' as the newline code. (Warning: Windows OS uses 'CR+LF' as default. You should specify the newline code when you save the list file.)
## (2) end with a blank line (Add a linefeed after the last sample. Otherwise, the last sample will be ignored.)
##【重要】配列名一覧ファイルでは
## (1) 改行コードは「LF」と指定すること （Windowsでは CR+LFなので，変更して保存すること）
## (2) 最後は空行とすること（最後のサンプル名の後にも改行を忘れずに入れること。そうしないと最後のサンプルがスキップされます）
#----------------------------------------

# Directories and files
## Suffix of the quality-trimmed sequence files / トリミング済配列ファイルの語尾
sfx_r1="_R1.qt.fq" # suffix of read 1 of the quality paired-end sequence
sfx_r2="_R2.qt.fq" # suffix of read 2 of the quality paired-end sequence

## the directory of quality-trimmed fastq files (specify the absolute path. do not include the final '/')  
## トリミング後のFastq配列データのあるディレクトリ（絶対パスで指定。最後のスラッシュ '/' は含めない）
DIR_SEQ="/home/ryohonda/sequence"
DIR_QT="${DIR_SEQ}/1.qt" # directory of fastp reports
DIR_ARG="${DIR_WORKING}/4.ARG_profile" # directory of ARG profile of each sample

## output directory 
## 結果出力するディレクトリ（絶対パスで指定。最後のスラッシュ '/' は含めない）
DIR_16S="${DIR_WORKING}/5.16S_taxonomy" # 16S
DIR_SUM="${DIR_WORKING}/8.summary" # summary of reads

## Name of database 検索するデータベース名
## location of the database データベースのあるディレクトリ
DIR_DB="/home/ryohonda/db/SILVA-138_k2"
## location of the python script
PY_SUM="${DIR_WORKING}/scripts/summary_reads_16S_ARG.py"

#====== Singularity (AppContainer) settings ==============================
# Singularity: ** specify the directory of other users if you need to refer from singularity.
# Singularity: **他のユーザのファイルをsingularityから参照する場合は次で指定 **必要ない場合は削除**
export SINGULARITY_BINDPATH="${DIR_DB},${DIR_SEQ}"
## Singularity: location of the container images of the singularity package
KRAKEN2="/usr/local/biotools/k/kraken2:2.1.2--pl5321h9f5acd7_3"
BRACKEN="/usr/local/biotools/b/bracken:2.8--py39hc16433a_0"
#==========================================================

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
# Execute from Singularity. Singularityから実行
# Check the file names in the command line below are correct. 下行のコマンドライン内のファイル名が正しいか確認。
## record the start of the script
DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
i=0; n=${#LIST[@]}
echo "[${DATE}] $JOB_NAME started. (${i}/${n})"
## create output directories if not existed
mkdir -p ${DIR_16S} ${DIR_SUM}
dir_rep="${DIR_16S}/reports"; mkdir -p ${dir_rep}; trap 'rm -R ${dir_rep}' 0


## Run for each file
for SAMPLE in "${LIST[@]}"; do
	## unzip fastq files 
	#dir_tmp="${DIR_QT}/tmp"; mkdir -p ${dir_tmp}; trap 'rm -R ${dir_tmp}' 0
	#gunzip -c ${DIR_QT}/${SAMPLE}_R1.qt.fq.gz > ${DIR_QT}/${SAMPLE}_R1.qt.fq
	#gunzip -c ${DIR_QT}/${SAMPLE}_R2.qt.fq.gz > ${DIR_QT}/${SAMPLE}_R2.qt.fq
	
	## kraken2 - for mpa-style report
	singularity exec ${KRAKEN2} \
	kraken2 --db ${DIR_DB} --threads ${threads} --use-mpa-style\
	 --report ${DIR_16S}/${SAMPLE}.mpa.tsv\
	 --paired ${DIR_QT}/${SAMPLE}${sfx_r1} ${DIR_QT}/${SAMPLE}${sfx_r2}\
	 > ${dir_rep}/${SAMPLE}.kraken2
	## kraken2 - default-style report
	singularity exec ${KRAKEN2} \
	kraken2 --db ${DIR_DB} --threads ${threads}\
	 --report ${dir_rep}/${SAMPLE}.report.kr2\
	 --paired ${DIR_QT}/${SAMPLE}${sfx_r1} ${DIR_QT}/${SAMPLE}${sfx_r2}\
	 > ${dir_rep}/${SAMPLE}.kraken2
	## delete unzipped fastq files
	#rm ${DIR_QT}/${SAMPLE}${sfx_r1} ${DIR_QT}/${SAMPLE}_R2.qt.fq

	## bracken - genus level
	singularity exec ${BRACKEN} \
	bracken -d ${DIR_DB} -i ${dir_rep}/${SAMPLE}.report.kr2\
	 -o ${DIR_16S}/${SAMPLE}_5.genus.tsv\
	 -r $rlen -l G
	## bracken - family level
	singularity exec ${BRACKEN} \
	bracken -d ${DIR_DB} -i ${dir_rep}/${SAMPLE}.report.kr2\
	 -o ${DIR_16S}/${SAMPLE}_4.family.tsv\
	 -r $rlen -l F
	## bracken - order level
	singularity exec ${BRACKEN} \
	bracken -d ${DIR_DB} -i ${dir_rep}/${SAMPLE}.report.kr2\
	 -o ${DIR_16S}/${SAMPLE}_3.order.tsv\
	 -r $rlen -l O
	## bracken - class level
	singularity exec ${BRACKEN} \
	bracken -d ${DIR_DB} -i ${dir_rep}/${SAMPLE}.report.kr2\
	 -o ${DIR_16S}/${SAMPLE}_2.class.tsv\
	 -r $rlen -l C
	## bracken - phylum level
	singularity exec ${BRACKEN} \
	bracken -d ${DIR_DB} -i ${dir_rep}/${SAMPLE}.report.kr2\
	 -o ${DIR_16S}/${SAMPLE}_1.phylum.tsv\
	 -r $rlen -l P
	## bracken - domain level
	singularity exec ${BRACKEN} \
	bracken -d ${DIR_DB} -i ${dir_rep}/${SAMPLE}.report.kr2\
	 -o ${DIR_16S}/${SAMPLE}_0.domain.tsv\
	 -r $rlen -l D
	
	## record the progress
	DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
	i=`expr $i + 1`
	echo "[${DATE}] ${SAMPLE} is done. (${i}/${n})"
done

## summarize reads data
python3 ${PY_SUM} ${DIR_QT} ${DIR_16S} ${DIR_ARG} ${DIR_SUM}

## record the completion of the script
DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
echo "[$DATE] $JOB_NAME completed. (${i}/${n})" 
exit 0
