#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -o $HOME/log/
#$ -e $HOME/log/
#$ -pe def_slot 8

# ** specify the directory of other users if you need to refer from singularity.
# **他のユーザの配列ファイルをsingularityから参照する場合は次で指定 **必要ない場合は削除**
export SINGULARITY_BINDPATH="/home/user/db,/home/user/JST-MIRAI"

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

## the directory of raw sequence files (specify the absolute path. do not include the final '/')  
## 配列データのあるディレクトリ（絶対パスで指定。最後のスラッシュ '/' は含めない）
DIR_SEQ="/home/user/JST-MIRAI"
DIR_QT="${DIR_SEQ}/1.trimmed" # directory of trimmed sequence and fastp reports
DIR_ARG="${DIR_SEQ}/5.ARG_profile" # directory of ARG profile of each sample
rlen=150 # read length

## output directory 
## 結果出力するディレクトリ（絶対パスで指定。最後のスラッシュ '/' は含めない）
DIR_16S="${DIR_WORKING}/6.16S" # 16S
DIR_SUM="${DIR_WORKING}/7.summary" # summary of reads

## Name of database 検索するデータベース名
## location of the database データベースのあるディレクトリ
DIR_DB="/home/user/db/SILVA-138_k2"

## location of the target package in singularity 
KRAKEN2="/usr/local/biotools/k/kraken2:2.1.2--pl5321h9f5acd7_3"
BRACKEN="/usr/local/biotools/b/bracken:2.8--py39hc16433a_0"

## location of the python script
PY_SUM="${DIR_WORKING}/script/summary_reads_16S_ARG.py"

####### Run program #############################################################
## Listing the raw sequence files.
if [ $LISTING -eq 2 ]; then 
	# in case of sequential files
	for ((i=0;i<=`expr $END - $START`;i++)); do
		n=`expr $i + $START`; LIST[$i]="${PREFIX}${n}"
	done
else
	# get the list from a file. 
	LIST=(`cat ${FILE_LIST}`)
fi

## Execute from Singularity. Singularityから実行
## Check the file names in the command line below are correct. 下行のコマンドライン内のファイル名が正しいか確認。
DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
i=0; n=${#LIST[@]}
echo "[${DATE}] $JOB_NAME started. (${i}/${n})"
## make output directories if not existed
dir_rep="${DIR_16S}/reports"; mkdir -p ${DIR_16S} ${dir_rep} ${DIR_SUM}

for SAMPLE in "${LIST[@]}"; do
	## unzip fastq files 
	#dir_tmp="${DIR_QT}/tmp"; mkdir -p ${dir_tmp}; trap 'rm -R ${dir_tmp}' 0
	#gunzip -c ${DIR_QT}/${SAMPLE}_R1.qt.fq.gz > ${DIR_QT}/${SAMPLE}_R1.qt.fq
	#gunzip -c ${DIR_QT}/${SAMPLE}_R2.qt.fq.gz > ${DIR_QT}/${SAMPLE}_R2.qt.fq
	
	## kraken2 - for mpa-style report
	singularity exec ${KRAKEN2} \
	kraken2 --db ${DIR_DB} --threads 8 --use-mpa-style\
	 --report ${DIR_16S}/${SAMPLE}.mpa.tsv\
	 --paired ${DIR_QT}/${SAMPLE}_R1.qt.fq ${DIR_QT}/${SAMPLE}_R2.qt.fq\
	 > ${dir_rep}/${SAMPLE}.kraken2
	## kraken2 - default-style report
	singularity exec ${KRAKEN2} \
	kraken2 --db ${DIR_DB} --threads 8\
	 --report ${dir_rep}/${SAMPLE}.report.kr2\
	 --paired ${DIR_QT}/${SAMPLE}_R1.qt.fq ${DIR_QT}/${SAMPLE}_R2.qt.fq\
	 > ${dir_rep}/${SAMPLE}.kraken2
	 
	## delete unzipped fastq files
	#rm ${DIR_QT}/${SAMPLE}_R1.qt.fq ${DIR_QT}/${SAMPLE}_R2.qt.fq

	## bracken - genus level
	singularity exec ${BRACKEN} \
	bracken -d ${DIR_DB} -i ${dir_rep}/${SAMPLE}.report.kr2\
	 -o ${DIR_16S}/${SAMPLE}_5.genus.tsv\
	 -w ${dir_rep}/${SAMPLE}.report_5.genus.bk.kr2\
	 -r $rlen -l G
	## bracken - family level
	singularity exec ${BRACKEN} \
	bracken -d ${DIR_DB} -i ${dir_rep}/${SAMPLE}.report.kr2\
	 -o ${DIR_16S}/${SAMPLE}_4.family.tsv\
	 -w ${dir_rep}/${SAMPLE}.report_4.family.bk.kr2\
	 -r $rlen -l F
	## bracken - order level
	singularity exec ${BRACKEN} \
	bracken -d ${DIR_DB} -i ${dir_rep}/${SAMPLE}.report.kr2\
	 -o ${DIR_16S}/${SAMPLE}_3.order.tsv\
	 -w ${dir_rep}/${SAMPLE}.report_3.order.bk.kr2\
	 -r $rlen -l O
	## bracken - class level
	singularity exec ${BRACKEN} \
	bracken -d ${DIR_DB} -i ${dir_rep}/${SAMPLE}.report.kr2\
	 -o ${DIR_16S}/${SAMPLE}_2.class.tsv\
	 -w ${dir_rep}/${SAMPLE}.report_2.class.bk.kr2\
	 -r $rlen -l C
	## bracken - phylum level
	singularity exec ${BRACKEN} \
	bracken -d ${DIR_DB} -i ${dir_rep}/${SAMPLE}.report.kr2\
	 -o ${DIR_16S}/${SAMPLE}_1.phylum.tsv\
	 -w ${dir_rep}/${SAMPLE}.report_1.phylum.bk.kr2\
	 -r $rlen -l P
	## bracken - domain level
	singularity exec ${BRACKEN} \
	bracken -d ${DIR_DB} -i ${dir_rep}/${SAMPLE}.report.kr2\
	 -o ${DIR_16S}/${SAMPLE}_0.domain.tsv\
	 -w ${dir_rep}/${SAMPLE}.report_0.domain.bk.kr2\
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
