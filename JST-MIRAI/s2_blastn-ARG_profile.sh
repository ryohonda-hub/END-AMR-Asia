#!/bin/sh
#==============================================================================
# blastn-ARG_profile, JST-MIRAI version / created by Ryo Honda (END-AMR-Asia)
#==============================================================================
# This shell script creates ARG profile data including:
#  - gene symbol and family, read counts, RPK, drug class and resistance mechanism
# This shell script recalls and requires: 
# 1. blastn [to blast the sequences on CARD database], 
# 2. count_blast_hits.pl [to count the reads of each ARG in the blast hits with >=90% of pident (% of identical positions) to >=25 bp of length], 
# 3. lookup_gene_info.py [to look up the gene information from the CARD catalog, and add in the read count data.]
#------------------------------------------------------------------------------

####### Parameter setting ######################################################
# Sample name prefix of sequence files, 配列ファイル名の最初についているサンプル名
SAMPLE="Prefix_"
# Numbering of samples to be processed 処理する配列ファイルの番号
START=165
END=165

# Your working directory and the directory of query sequence files  (specify the absolute path. do not include the final '/')  
# 作業ディレクトリと配列データのあるディレクトリ（絶対パスで指定。最後のスラッシュ '/' は含めない）
DIR_WORKING="/home/user/JST-MIRAI"
DIR_QUERY="${DIR_WORKING}/2.fasta"

# Output directory  (specify the absolute path. do not include the final '/') 
# 結果出力するディレクトリ（絶対パスで指定。最後のスラッシュ '/' は含めない）
DIR_BLAST="${DIR_WORKING}/3.blast_CARD"
DIR_COUNT="${DIR_WORKING}/4.count"
DIR_ARG="${DIR_WORKING}/5.ARG_profile"

# Name of database to be specified in blast command 検索するデータベース名 (blastn で指定する値）
# location of the database and ARO gene catalog データベースのあるディレクトリと，遺伝子情報カタログ(aro_index)
DB="CARD-3.0.7"
DIR_DB="/home/user/db/CARD-3.0.7"
DB_CAT="aro_index.tsv"

# location of the scripts for counting hits and matching gene information
PL_COUNT="/home/user/scripts/count_blast_hits.pl"
PY_LKUP="/home/user/scripts/make_ARG_profile.py"

####### Run program #############################################################
# Execute from Singularity. Singularityから実行
# Check the file names in the command line below are correct. 下行のコマンドライン内のファイル名が正しいか確認。

DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
echo "[${DATE}] $JOB_NAME started. ${SAMPLE}${START}-${END}"
# make output directories if not existed
mkdir -p ${DIR_BLAST} "${DIR_BLAST}/filtered" ${DIR_COUNT} ${DIR_ARG}

for ((i=START; i<=END; i++))
do
	# blastn
	blastn -query ${DIR_QUERY}/${SAMPLE}${i}_R12.fa -db ${DIR_DB}/${DB}\
	-out ${DIR_BLAST}/${SAMPLE}${i}.blast.txt -num_threads 24\
	-max_target_seqs 1 -evalue 1E-5\
	-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
	
	# exclude the blast results with pident <90% or length <25bp
	awk '{if($3 >= 90 && $4 >= 25) print $0}' ${DIR_BLAST}/${SAMPLE}${i}.blast.txt > ${DIR_BLAST}/filtered/${SAMPLE}${i}.blast.filtered.txt
	# count gene hits
	perl ${PL_COUNT} ${DIR_BLAST}/filtered/${SAMPLE}${i}.blast.filtered.txt > ${DIR_COUNT}/${SAMPLE}${i}.count.tsv
	# look up gene information in the CARD catalog
	python ${PY_LKUP} ${DIR_DB}/${DB_CAT} ${DIR_COUNT}/${SAMPLE}${i}.count.tsv ${DIR_ARG}
	
	# record the progress
	DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
	echo "[${DATE}] ${SAMPLE}${i} is done."
done

# record the completion of the script
DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
echo "[$DATE] $JOB_NAME completed."

exit 0
