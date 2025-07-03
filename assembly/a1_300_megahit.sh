#!/bin/bash
JOB_NAME="${SLURM_JOB_NAME}" # any job name for your reference
#==============================================================================
# Assembly by Megahit, version 1.0 , Last updated: 2025-04-28
# created by Mardalisa and Ryo Honda
# The shell script for NIG supercomputer using AppContainer (Singularity)
#==============================================================================
# Threads recommended: 8-16 threads
# Memory recommended: 8GB per thread
#------------------------------------------------------------------------------

threads=8 # CPU threads / スレッド数
####### Parameter setting ######################################################
# Directories and files (Specify a directory as a full path without the last '/').
## Your working directory 自分の作業ディレクトリ
DIR_WORKING="/home/ryohonda/GlobalAMR/Assembly"
## The directory of the input sequence files / 配列ファイルのあるディレクトリ
DIR_SEQ="/home/ryohonda/GlobalAMR"
DIR_QT="${DIR_SEQ}/1.quality/Asia"
## Suffix of the input fastq files / アセンブリするFastq配列ファイルの語尾
sfx_r1="_R1.qt.fq"
sfx_r2="_R2.qt.fq"

## Output directory for contig files  
## 出力ディレクトリ（絶対パスで指定。最後のスラッシュ '/' は含めない）
DIR_CTG="${DIR_WORKING}/2.megahit/Asia300"

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
FILE_LIST="${DIR_WORKING}/sh/sralist.txt"
## [IMPORTANT] In the sample list file,
## (1) specify 'LF' as the newline code. (Warning: Windows OS uses 'CR+LF' as default. You should specify the newline code when you save the list file.)
## (2) end with a blank line (Add a linefeed after the last sample. Otherwise, the last sample will be ignored.)
##【重要】配列名一覧ファイルでは
## (1) 改行コードは「LF」と指定すること （Windowsでは CR+LFなので，変更して保存すること）
## (2) 最後は空行とすること（最後のサンプル名の後にも改行を忘れずに入れること。そうしないと最後のサンプルがスキップされます）
#----------------------------------------

# Assembly options
MIN_CTG=300
## location of the python script
PY_CHK="${DIR_WORKING}/scripts/chfa.py"
#====== Singularity (AppContainer) settings ==============================
# Singularity: ** specify the directory of other users if you need to refer from singularity.
# Singularity: **他のユーザのファイルをsingularityから参照する場合は次で指定 **必要ない場合は削除**
export SINGULARITY_BINDPATH="${DIR_QT}"
## Singularity: location of the container image of the singularity package
# location of the container image for singularity 
MEGAHIT="/usr/local/biotools/m/megahit:1.2.9--h8b12597_0"
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
DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
i=0; n=${#LIST[@]}
echo "[${DATE}] $JOB_NAME started. (${i}/${n})"
## create output directories if not existed
mkdir -p ${DIR_CTG}

for SAMPLE in "${LIST[@]}"; do
	# record the progress
	DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
	i=`expr $i + 1`
	echo "[${DATE}] ${SAMPLE} is started. (${i}/${n})"
	
	# run from the AppContainer (Singularity)
	singularity exec ${MEGAHIT} \
	megahit -1 ${DIR_QT}/${SAMPLE}${sfx_r1} -2 ${DIR_QT}/${SAMPLE}${sfx_r2}\
	 -o ${DIR_CTG}/${SAMPLE} --out-prefix ${SAMPLE}\
	 --min-contig-len ${MIN_CTG}\
	 -t ${threads} 

	# check contig length
	python3 ${PY_CHK} ${DIR_CTG}/${SAMPLE}/${SAMPLE}.contigs.fa > ${DIR_CTG}/${SAMPLE}.ctg.summary.txt
	 
	# record the progress
	DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
	echo "[${DATE}] ${SAMPLE} is done. (${i}/${n})"
done
# record the completion of the script
DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
echo "[$DATE] $JOB_NAME completed. (${i}/${n})"
exit 0
