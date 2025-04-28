#!/bin/bash
JOB_NAME="${SLURM_JOB_NAME}" # any job name for your reference
#==============================================================================
# create ARG_profile from blast hits data / created by Ryo Honda, Last updated: 2025-03-15
# The shell script for NIG supercomputer with the Singularity package of BioContainer
#==============================================================================
# This shell script creates ARG profile data from blast hit data including:
#  - gene symbol and family, read counts, RPK, drug class and resistance mechanism
# This shell script recalls and requires: 
# - makeARGprof.py [to count the reads of each ARG in the blast hits with pident (% of identical positions) or matching length below the cutoff thresholds, and to look up the gene information from the CARD catalog, and add in the read count data.]
#
# * blast output format should be:
# 	-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
#------------------------------------------------------------------------------

####### Parameter setting ######################################################
# Your working directory 自分の作業ディレクトリ
DIR_WORKING="/home/ryohonda/GlobalAMR"

# Sample list
## Choose listing method of sequence data files (choose 0 or 1)
## 配列ファイルリストの指定方法 (0か1を選択）
LISTING="0" 
#----------------------------------------
# 0: to use all the raw files in the directory specified as DIR_BLAST below.
#	 DIR_BLASTで指定したディレクトリーのすべてのファイルを使用する。
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
## Suffix of blast output files / Blast結果ファイルの語尾
sfx=".blast.txt" # suffix of the query fasta files

## The directory of blast output files
## blast結果データのあるディレクトリ（絶対パスで指定。最後のスラッシュ '/' は含めない）
DIR_BLAST="${DIR_WORKING}/3.blast_CARD"

## Output directory  (specify the absolute path. do not include the final '/') 
## 結果出力するディレクトリ（絶対パスで指定。最後のスラッシュ '/' は含めない）
DIR_PROF="${DIR_WORKING}/4.ARG_profile"

## Name of database to be specified in blast command 検索するデータベース名 (blastn で指定する値）
## location of the database and ARO gene catalog データベースのあるディレクトリと，遺伝子情報カタログ(aro_index)
DB="CARD-3.2.6"
DIR_DB="/home/ryohonda/db/CARD-3.2.6"
DB_CAT="aro_index.tsv"
## location of the scripts for counting hits and matching gene information
PY_LKUP="${DIR_WORKING}/scripts/makeARGprof.py"

#====== Singularity (AppContainer) settings ==============================
# Singularity: ** specify the directory of other users if you need to refer from singularity.
# Singularity: **他のユーザのファイルをsingularityから参照する場合は次で指定 **必要ない場合は削除**
export SINGULARITY_BINDPATH="${DIR_DB},${DIR_BLAST}"
## Singularity: location of the container image of the singularity package
BLAST="/usr/local/biotools/b/blast:2.9.0--pl526he19e7b1_7"
#==========================================================

####### Run program #############################################################
## Listing the raw sequence files.
if [ $LISTING -eq 0 ]; then 
	# Case 0: get the sample list from the fastq file list of the directory 
	LIST=($(ls ${DIR_BLAST}/*${sfx} | sed -E "s/${sfx}$//" | xargs -n 1 basename | sort -u))
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
	# Create ARG Profile
	## exclude the blast results with low pident or short matching length
	## count gene hits
	## look up gene information in the CARD catalog
	python3 ${PY_LKUP} ${DIR_DB}/${DB_CAT} ${DIR_BLAST}/${SAMPLE}${sfx}  ${DIR_PROF}
	
	## record the progress
	DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
	i=`expr $i + 1`
	echo "[${DATE}] ${SAMPLE} is done. (${i}/${n})"
done
## record the completion of the script
DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
echo "[$DATE] $JOB_NAME completed. (${i}/${n})"
exit 0