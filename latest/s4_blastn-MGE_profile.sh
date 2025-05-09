#!/bin/bash
JOB_NAME="${SLURM_JOB_NAME}" # any job name for your reference
#==============================================================================
# blastn-MBE_profile, version 2.1 / created by Ryo Honda, Last updated: 2025-03-15
# The shell script for NIG supercomputer with the Singularity package of BioContainer
#==============================================================================
# This shell script creates MGE profile data including:
#  - gene symbol and family, read counts, RPK, drug class and resistance mechanism
# This shell script recalls and requires: 
# 1. blastn [to identify MGEs using the MGE database], 
# 2. makeMGEprof.py [to count the reads of each MGE in the blast hits with pident (% of identical positions) or matching length below the cutoff thresholds, and to look up the gene information from the MGEDB catalog, and add in the read count data.]
#------------------------------------------------------------------------------

threads=24	# CPU threads / スレッド数
####### Parameter setting ######################################################
# Your working directory 自分の作業ディレクトリ
DIR_WORKING="/home/ryohonda/GlobalAMR"

# Sample list
## Choose listing method of sequence data files (choose 0 or 1)
## 配列ファイルリストの指定方法 (0か1を選択）
LISTING="1" 
#----------------------------------------
# 0: to use all the raw files in the directory specified as DIR_QUERY below.
#	 DIR_QUERYで指定したディレクトリーのすべてのファイルを使用する。
#----------------------------------------
# 1: to input the file list from a text file (1 file in 1 line).
#    配列リストをファイルから読み込む（1行1ファイルで記述）
#
## For Case 1: Specify the text file of sequence list and newline code
## 1の場合: 配列名一覧のファイルと改行コードを指定
FILE_LIST="${DIR_WORKING}/sralist.txt"
## [IMPORTANT] In the sample list file,
## (1) specify 'LF' as the newline code. (Warning: Windows OS uses 'CR+LF' as default. You should specify the newline code when you save the list file.)
## (2) end with a blank line (Add a linefeed after the last sample. Otherwise, the last sample will be ignored.)
##【重要】配列名一覧ファイルでは
## (1) 改行コードは「LF」と指定すること （Windowsでは CR+LFなので，変更して保存すること）
## (2) 最後は空行とすること（最後のサンプル名の後にも改行を忘れずに入れること。そうしないと最後のサンプルがスキップされます）
#----------------------------------------

# Directories and files
## Your working directory and the directory of query sequence files  (specify the absolute path. do not include the final '/')  
## Suffix of query sequence files / QueryとなるFasta配列ファイルの語尾
sfx_fa="_R12.fa" # suffix of query fasta files

## 配列データのあるディレクトリ（絶対パスで指定。最後のスラッシュ '/' は含めない）
DIR_SEQ="/home/ryohonda/sequence"
DIR_QUERY="${DIR_SEQ}/2.fasta"
## Output directory  (specify the absolute path. do not include the final '/') 
## 結果出力するディレクトリ（絶対パスで指定。最後のスラッシュ '/' は含めない）
DIR_BLAST="${DIR_WORKING}/6.blast_MGEDB"
DIR_PROF="${DIR_WORKING}/7.MGE_profile"

## Name of database to be specified in blast command 検索するデータベース名 (blastn で指定する値）
## location of the database and MGE gene catalog データベースのあるディレクトリと，遺伝子情報カタログ(MGE_tax_table_trim.txt)
DB="MGEDB"
DIR_DB="/home/ryohonda/db/MGEDB"
DB_CAT="MGE_tax_table_trim.txt"
## location of the scripts for counting hits and matching gene information
PY_LKUP="${DIR_WORKING}/scripts/makeMGEprof.py"

#====== Singularity  (AppContainer) settings ==============================
# Singularity: ** specify the directory of other users if you need to refer from singularity.
# Singularity: **他のユーザのファイルをsingularityから参照する場合は次で指定 **必要ない場合は削除**
export SINGULARITY_BINDPATH="${DIR_DB},${DIR_SEQ}"
## Singularity: location of the container image of the singularity package
BLAST="/usr/local/biotools/b/blast:2.9.0--pl526he19e7b1_7"
#==========================================================

####### Run program #############################################################
## Listing the raw sequence files.
if [ $LISTING -eq 0 ]; then 
	# Case 0: get the sample list from the fastq file list of the directory 
	LIST=($(ls ${DIR_QUERY}/*${sfx_fa} | sed -E "s/${sfx_fa}$//" | xargs -n 1 basename | sort -u))
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
