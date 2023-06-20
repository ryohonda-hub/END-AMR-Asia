#! /bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -o $HOME/log/
#$ -e $HOME/log/
#$ -pe def_slot 16

# ** specify the directory of other users if you need to refer.
# **他のユーザのファイルを参照する場合は次で指定 **必要ない場合は削除**
export SINGULARITY_BINDPATH="/home/user/JST-MIRAI"

####### Parameter setting ######################################################
## NCBIからダウンロードするfastqファイルの保存先
DIR_OUT="/home/user/JST-MIRAI/0.raw"

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

# location of the container image for singularity 
SRATOOL="/usr/local/biotools/s/sra-tools:3.0.5--h9f5acd7_1"

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

## Execute from Singularity. Singularityから実行
DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
i=0; n=${#LIST[@]}
echo "[${DATE}] $JOB_NAME started. (${i}/${n})"

for DRR in "${LIST[@]}"; do
	DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
	echo "[$DATE] $DRR started."
	if[[-e "${DIR_OUT}/${DRR}_1.fastq"]]; then
		echo "${DRR} was not downloaded because it already exists in ${DIR_OUT}."
	else
		singularity exec ${SRATOOL} \
		fasterq-dump ${DRR} -e 16 -t ${DIR_OUT} --outdir ${DIR_OUT}
		DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
		echo "[$DATE] $DRR done."
	fi
	## record the progress
	DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
	i=`expr $i + 1`
	echo "[${DATE}] ${SAMPLE} is done. (${i}/${n})"
done
## record the completion of the script
DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
echo "[$DATE] ${JOB_NAME} completed. ${i} files were downloaded in ${DIR_OUT}."

