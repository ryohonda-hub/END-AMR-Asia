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
## Directory of input contig files / 解析するコンティグファイルのディレクトリ
DIR_CTG="${DIR_WORKING}/2.megahit/Asia" 
## Suffix of the contig files / アセンブリしたコンティグファイルの語尾
sfx_fa=".contigs.fa" # suffix of the contig files

## Output directory 
## 結果出力するディレクトリ（絶対パスで指定。最後のスラッシュ '/' は含めない）
DIR_ARG="${DIR_WORKING}/3.ARGI/Asia" # RGI results of contigs
DIR_PROF="${DIR_WORKING}/4.ARGI_profile/Asia" # ARG profile of contigs
DIR_16S="${DIR_WORKING}/5.16S_taxonomy/Asia" # 16S taxonomy of contigs
DIR_SUM="${DIR_WORKING}/6.summary" # ARG-host table

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
DIR_SILVA="/home/ryohonda/db/SILVA-138_k2"
DIR_TAXDMP="/home/ryohonda/db/taxdump"
DB_TAXDMP="names.dmp"
DIR_GTDB="/home/ryohonda/db/kraken2/GTDBv226"
DB_GTDB="inspect.txt"
## location of the python script
PY_PROF="${DIR_WORKING}/scripts/ARGIprof.py"
PY_HOST="${DIR_WORKING}/scripts/ARGIhost.py"

#====== Singularity (AppContainer) settings ==============================
# Singularity: ** specify the directory of other users if you need to refer from singularity.
# Singularity: **他のユーザのファイルをsingularityから参照する場合は次で指定 **必要ない場合は削除**
export SINGULARITY_BINDPATH="${DIR_CARD},${DIR_SILVA},${DIR_TAXDMP},${DIR_GTDB}"
## Singularity: location of the container images of the singularity package
RGI="/usr/local/biotools/r/rgi:6.0.3--pyha8f3691_1"
KRAKEN2="/usr/local/biotools/k/kraken2:2.1.2--pl5321h9f5acd7_3"
BRACKEN="/usr/local/biotools/b/bracken:2.8--py39hc16433a_0"
#==========================================================

####### Run program #############################################################
## Listing the sample names.
if [ $LISTING -eq 0 ]; then 
	# Case 0: get the sample list from the directory list 
	LIST=($(find "${DIR_CTG}" -mindepth 1 -maxdepth 1 -type d -exec basename {} \;))
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
mkdir -p ${DIR_ARG} ${DIR_PROF} ${DIR_16S} ${DIR_SUM}
dir_rep="${DIR_16S}/reports"; mkdir -p ${dir_rep}; trap 'rm -R ${dir_rep}' 0

# Load CARD database
singularity exec ${RGI} \
rgi load --card_json ${DIR_CARD}/${DB_CARD} --local

## Run for each file
for SAMPLE in "${LIST[@]}"; do
	## record the progress
	DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
	i=`expr $i + 1`
	echo "[${DATE}] ${SAMPLE}: Identifying ARG by RGI started. (${i}/${n})"

	# Identification of ARGs on contigs
	singularity exec ${RGI} \
	rgi main\
     -i ${DIR_CTG}/${SAMPLE}/${SAMPLE}${sfx_fa}\
     -o ${DIR_ARG}/${SAMPLE}\
     -t contig -a DIAMOND -n ${threads}\
     --local --clean --low_quality --include_loose
	
	## record the progress
	DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
	echo "[${DATE}] ${SAMPLE}: Identifying ARG by RGI is done."
	echo "[${DATE}] ${SAMPLE}: Creating the ARG profile table."

	## Create the ARG profile table
	python3 ${PY_PROF} ${DIR_CARD}/${DB_CARD} ${DIR_ARG}/${SAMPLE}.ARG_profile.tsv ${DIR_PROF}
	
	## record the progress
	DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
	echo "[${DATE}] ${SAMPLE}: ARG profiling is done."
	echo "[${DATE}] ${SAMPLE}: Phylogenetic analysis by kraken2 started."

	# Phylegentic analysis by kraken2+bracken
	## kraken2 - for mpa-style report
	singularity exec ${KRAKEN2} \
	kraken2 --db ${DIR_SILVA} --threads ${threads} --use-mpa-style\
	 --report ${DIR_16S}/${SAMPLE}.mpa.tsv\
	 ${DIR_CTG}/${SAMPLE}/${SAMPLE}${sfx_fa}\
	 > ${DIR_16S}/${SAMPLE}.mpa.kraken2.tsv
	## kraken2 - default-style report
	singularity exec ${KRAKEN2} \
	kraken2 --db ${DIR_SILVA} --threads ${threads}\
	 --report ${dir_rep}/${SAMPLE}.report.kr2\
	 ${DIR_CTG}/${SAMPLE}/${SAMPLE}${sfx_fa}\
	 > ${DIR_16S}/${SAMPLE}.kraken2.tsv

	## bracken - specie level
	singularity exec ${BRACKEN} \
	bracken -d ${DIR_SILVA} -i ${dir_rep}/${SAMPLE}.report.kr2\
	 -o ${DIR_16S}/${SAMPLE}_6.specie.tsv\
	 -l S
	## bracken - genus level
	singularity exec ${BRACKEN} \
	bracken -d ${DIR_SILVA} -i ${dir_rep}/${SAMPLE}.report.kr2\
	 -o ${DIR_16S}/${SAMPLE}_5.genus.tsv\
	 -l G
	## bracken - family level
	singularity exec ${BRACKEN} \
	bracken -d ${DIR_SILVA} -i ${dir_rep}/${SAMPLE}.report.kr2\
	 -o ${DIR_16S}/${SAMPLE}_4.family.tsv\
	 -l F
	## bracken - order level
	singularity exec ${BRACKEN} \
	bracken -d ${DIR_SILVA} -i ${dir_rep}/${SAMPLE}.report.kr2\
	 -o ${DIR_16S}/${SAMPLE}_3.order.tsv\
	 -l O
	## bracken - class level
	singularity exec ${BRACKEN} \
	bracken -d ${DIR_SILVA} -i ${dir_rep}/${SAMPLE}.report.kr2\
	 -o ${DIR_16S}/${SAMPLE}_2.class.tsv\
	 -l C
	## bracken - phylum level
	singularity exec ${BRACKEN} \
	bracken -d ${DIR_SILVA} -i ${dir_rep}/${SAMPLE}.report.kr2\
	 -o ${DIR_16S}/${SAMPLE}_1.phylum.tsv\
	 -l P
	## bracken - domain level
	singularity exec ${BRACKEN} \
	bracken -d ${DIR_SILVA} -i ${dir_rep}/${SAMPLE}.report.kr2\
	 -o ${DIR_16S}/${SAMPLE}_0.domain.tsv\
	 -l D
	
	## record the progress
	DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
	echo "[${DATE}] ${SAMPLE}: Phylogenetic analysis by kraken2 is done. (${i}/${n})"
	
	# Create ARG-host table
	python3 ${PY_HOST} ${DIR_GTDB}/${DB_GTDB} ${DIR_ARG}/${SAMPLE}_ARGI.txt ${DIR_16S}/${SAMPLE}.kraken2.tsv ${DIR_SUM}
done
# Unload CARD database
singularity exec ${RGI} \
rgi clean --local

## record the completion of the script
DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
echo "[$DATE] $JOB_NAME completed. (${i}/${n})" 
exit 0
