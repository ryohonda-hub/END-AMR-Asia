#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -o $HOME/log/
#$ -e $HOME/log/
#$ -pe def_slot 24

#==============================================================================
# blastn-ARG_profile, JST-MIRAI+global version / created by Ryo Honda, 2023-06-21
# The shell script for NIG supercomputer (Grid Engine) 
# with Singularity package of BioContainer
#==============================================================================
# This shell script creates ARG profile data including:
#  - gene symbol and family, read counts, RPK, drug class and resistance mechanism
# This shell script recalls and requires: 
# 1. blastn [to blast the sequences on CARD database], 
# 2. count_blast_hits.pl [to count the reads of each ARG in the blast hits with >=90% of pident (% of identical positions) to >=25 bp of length], 
# 3. lookup_gene_info.py [to look up the gene information from the CARD catalog, and add in the read count data.]
#------------------------------------------------------------------------------

# ** specify the directory of other users if you need to refer from singularity.
# **他のユーザのファイルをsingularityから参照する場合は次で指定 **必要ない場合は削除**
export SINGULARITY_BINDPATH="/home/user/db"
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

# Your working directory and the directory of query sequence files  (specify the absolute path. do not include the final '/')  
# 配列データのあるディレクトリ（絶対パスで指定。最後のスラッシュ '/' は含めない）
DIR_SEQ="/home/user/JST-MIRAI"
DIR_QUERY="${DIR_SEQ}/2.fasta"

# Output directory  (specify the absolute path. do not include the final '/') 
# 結果出力するディレクトリ（絶対パスで指定。最後のスラッシュ '/' は含めない）
DIR_BLAST="${DIR_WORKING}/3.blast_CARD"
DIR_COUNT="${DIR_WORKING}/4.count"
DIR_ARG="${DIR_WORKING}/5.ARG_profile"

# Name of database to be specified in blast command 検索するデータベース名 (blastn で指定する値）
# location of the database and ARO gene catalog データベースのあるディレクトリと，遺伝子情報カタログ(aro_index)
DB="CARD-3.2.6"
DIR_DB="/home/user/db/CARD-3.2.6"
DB_CAT="aro_index.tsv"

# location of the scripts for counting hits and matching gene information
PL_COUNT="/home/user/JST-MIRAI/scripts/count_blast_hits.pl"
PY_LKUP="/home/user/JST-MIRAI/scripts/make_ARG_profile.py"
# location of the container image for singularity 
BLAST="/usr/local/biotools/b/blast:2.9.0--pl526he19e7b1_7"

####### Run program #############################################################
## Listing the raw sequence files.
if [ $LISTING -eq 1 ]; then 
	# in case of sequential files
	for ((i=0;i<=`expr $END - $START`;i++)); do
		n=`expr $i + $START`; LIST[$i]="${PREFIX}${n}"
	done
else
	# get the list from a file. 
	LIST=(`cat ${FILE_LIST}`)
fi

# Execute from Singularity. Singularityから実行
# Check the file names in the command line below are correct. 下行のコマンドライン内のファイル名が正しいか確認。

DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
i=0; n=${#LIST[@]}
echo "[${DATE}] $JOB_NAME started. (${i}/${n})"
# make output directories if not existed
mkdir -p ${DIR_BLAST} "${DIR_BLAST}/filtered" ${DIR_COUNT} ${DIR_ARG}

for SAMPLE in "${LIST[@]}"; do
	# blastn
	singularity exec ${BLAST} \
	blastn -query ${DIR_QUERY}/${SAMPLE}_R12.fa -db ${DIR_DB}/${DB}\
	-out ${DIR_BLAST}/${SAMPLE}.blast.txt -num_threads 24\
	-max_target_seqs 1 -evalue 1E-5\
	-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
	
	# exclude the blast results with pident <90% or length <25bp
	awk '{if($3 >= 90 && $4 >= 25) print $0}' ${DIR_BLAST}/${SAMPLE}.blast.txt > ${DIR_BLAST}/filtered/${SAMPLE}.blast.filtered.txt
	# count gene hits
	perl ${PL_COUNT} ${DIR_BLAST}/filtered/${SAMPLE}.blast.filtered.txt > ${DIR_COUNT}/${SAMPLE}.count.tsv
	# look up gene information in the CARD catalog
	python3 ${PY_LKUP} ${DIR_DB}/${DB_CAT} ${DIR_COUNT}/${SAMPLE}.count.tsv ${DIR_ARG}
	
	# record the progress
	DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
	i=`expr $i + 1`
	echo "[${DATE}] ${SAMPLE} is done. (${i}/${n})"
done

# record the completion of the script
DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
echo "[$DATE] $JOB_NAME completed. (${i}/${n})"

exit 0
