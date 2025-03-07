#!/bin/bash
#==============================================================================
# blastn-MBE_profile, version 2.1 / created by Ryo Honda, 2025-03-07
# The shell script for NIG supercomputer
# with Singularity package of BioContainer
#==============================================================================
# This shell script creates MGE profile data including:
#  - gene symbol and family, read counts, RPK, drug class and resistance mechanism
# This shell script recalls and requires: 
# 1. blastn [to blast the sequences on MGE database], 
# 2. make_MGE_profile2.py [to count the reads of each MGE in the blast hits with pident (% of identical positions) or matching length below the cutoff thresholds, and to look up the gene information from the MGEDB catalog, and add in the read count data.]
#------------------------------------------------------------------------------

threads=24	# CPU threads / スレッド数
####### Parameter setting ######################################################
# Your working directory 自分の作業ディレクトリ
DIR_WORKING="/home/ryohonda/GlobalAMR"

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

# Directories and files
## Your working directory and the directory of query sequence files  (specify the absolute path. do not include the final '/')  
## Suffix of raw sequence files / Raw配列ファイルの語尾
sfx_fa="_R12.fa" # suffix of read 1 of the paired-end sequence
## 配列データのあるディレクトリ（絶対パスで指定。最後のスラッシュ '/' は含めない）
DIR_SEQ="/home/ryohonda/sequence"
DIR_QUERY="${DIR_SEQ}/2.fasta"
## Output directory  (specify the absolute path. do not include the final '/') 
## 結果出力するディレクトリ（絶対パスで指定。最後のスラッシュ '/' は含めない）
DIR_BLAST="${DIR_WORKING}/7.blast_MGEDB"
DIR_PROF="${DIR_WORKING}/9.MGE_profile"

## Name of database to be specified in blast command 検索するデータベース名 (blastn で指定する値）
## location of the database and MGE gene catalog データベースのあるディレクトリと，遺伝子情報カタログ(MGE_tax_table_trim.txt)
DB="MGEDB"
DIR_DB="/home/ryohonda/db/MGEDB"
DB_CAT="MGE_tax_table_trim.txt"
## location of the scripts for counting hits and matching gene information
PY_LKUP="${DIR_WORKING}/scripts/make_MGE_prof2.py"

#====== Singularity settings ==============================
# Singularity: ** specify the directory of other users if you need to refer from singularity.
# Singularity: **他のユーザのファイルをsingularityから参照する場合は次で指定 **必要ない場合は削除**
export SINGULARITY_BINDPATH="${DIR_DB},${DIR_SEQ}"
## Singularity: location of the container image for singularity 
BLAST="/usr/local/biotools/b/blast:2.9.0--pl526he19e7b1_7"
#==========================================================

####### Run program #############################################################
## Listing the raw sequence files.
if [ $LISTING -eq 0 ]; then 
	# Case 0: get the sample list from the fastq file list of the directory 
	LIST=(`ls ${DIR_QUERY}/*${sfx_fa} | sed -E "s/${sfx_fa}$//" | sort -u`)
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
mkdir -p ${DIR_BLAST} ${DIR_PROF}

## Run for each file
for SAMPLE in "${LIST[@]}"; do
	# blastn
	singularity exec ${BLAST} \
	blastn -query ${DIR_QUERY}/${SAMPLE}${sfx_fa} -db ${DIR_DB}/${DB}\
	-out ${DIR_BLAST}/${SAMPLE}.blast.txt -num_threads ${threads}\
	-max_target_seqs 1 -evalue 1E-5\
	-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
	# Profiling
	## exclude the blast results with low pident or short matching length
	## count gene hits
	## look up gene information in the CARD catalog
	python3 ${PY_LKUP} ${DIR_DB}/${DB_CAT} ${DIR_BLAST}/${SAMPLE}.blast.txt ${DIR_PROF}
	
	## record the progress
	DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
	i=`expr $i + 1`
	echo "[${DATE}] ${SAMPLE} is done. (${i}/${n})"
done
## record the completion of the script
DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
echo "[$DATE] $JOB_NAME completed. (${i}/${n})"
exit 0
