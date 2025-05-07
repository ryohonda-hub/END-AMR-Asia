#!/bin/bash
JOB_NAME="${SLURM_JOB_NAME}" # any job name for your reference
#==============================================================================
# ARG host analysis from contigs by megahit, Last updated: 2025-05-01
# created by Ryo Honda
# The shell script for NIG supercomputer using AppContainer (Singularity)
#==============================================================================
# This shell script creates ARG-host profile tables.
#  - 16S-based taxonomy composition tables of each sample
#  - a table on the number of raw reads, quality reads, total 16S reads, and total ARG reads
# # This shell script recalls and requires: 
# 1. RGI and CARD database [for identifying ARG on contigs]
# 2. Kraken2 and Bracken2 [for taxonomy analysis using the SILVA database], 
# 3. summary_reads_16S_ARG.py [to creates a summary table of sequence reads]
#------------------------------------------------------------------------------

threads=8	# CPU threads / スレッド数
####### Parameter setting ######################################################
# Directories and files (Specify a directory as a full path without the last '/').
## Your working directory 自分の作業ディレクトリ
DIR_WORKING="/home/ryohonda/GlobalAMR/Assembly"

## Input and output directories
DIR_ARG="${DIR_WORKING}/3r.ARGI/Asia" # RGI results of contigs
DIR_PROF="${DIR_WORKING}/4.ARGI_profile/Asia" # ARG profile of contigs

## Suffix of rgit output files / Blast結果ファイルの語尾
sfx="_ARGI.txt" # suffix of the query fasta files

# Sample list
## Choose listing method of sequence data files (choose 0 or 1)
## 配列ファイルリストの指定方法 (0か1を選択）
LISTING="0" 
#----------------------------------------
# 0: to use all the sequence files in the directory specified as DIR_CTW above.
#	 DIR_CTGで指定したディレクトリーのすべてのサンプルを使用する。
#----------------------------------------
# 1: to input the file list from a text file (1 file in 1 line).
#    配列リストをファイルから読み込む（1行1ファイルで記述）
#
## For Case 1: Specify the text file of sequence list and newline code
## 1の場合: 配列名一覧のファイルと改行コードを指定
FILE_LIST="${DIR_WORKING}/sh/sralist.txt"
## [IMPORTANT] In the sample list file,
## (1) specify 'LF' as the newline code. (Warning: Windows OS uses 'CR+LF' as default. You should specify the newline code when you save the list file.)
## (2) end with a blank line (Add a linefeed after the last sample. Otherwise, the last sample will be ignored.)
##【重要】配列名一覧ファイルでは
## (1) 改行コードは「LF」と指定すること （Windowsでは CR+LFなので，変更して保存すること）
## (2) 最後は空行とすること（最後のサンプル名の後にも改行を忘れずに入れること。そうしないと最後のサンプルがスキップされます）
#----------------------------------------

## Name of database 検索するデータベース名
## location of the database データベースのあるディレクトリ
DIR_CARD="/home/ryohonda/db/CARD-3.2.6"
DB_CARD="card.json"	# specify CARD JSON file for RGI
## location of the python script
PY_PROF="${DIR_WORKING}/scripts/ARGIprof.py"

####### Run program #############################################################
## Listing the sample names.
if [ $LISTING -eq 0 ]; then 
	# Case 0: get the sample list from the directory list 
	LIST=($(ls ${DIR_ARG}/*${sfx} | sed -E "s/${sfx}$//" | xargs -n 1 basename | sort -u))
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
mkdir -p ${DIR_PROF} 

## Run for each file
for SAMPLE in "${LIST[@]}"; do
	## Create the ARG profile table
	python3 ${PY_PROF} ${DIR_CARD}/${DB_CARD} ${DIR_ARG}/${SAMPLE}${sfx} ${DIR_PROF}
	
	## record the progress
	DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
	i=`expr $i + 1`
	echo "[${DATE}] ${SAMPLE}: ARG profiling is done. (${i}/${n})"
done

## record the completion of the script
DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
echo "[$DATE] $JOB_NAME completed. (${i}/${n})" 
exit 0
