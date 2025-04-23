#!/bin/bash
JOB_NAME="${SLURM_JOB_NAME}" # any job name for your reference
#==============================================================================
# summarize sequence reads information (ARG, MGE and 16S) / created by Ryo Honda, Last updated: 2025-03-15
#==============================================================================
# This shell script creates a summary table of sequence reads of each sample including:
#  - the number of raw reads, quality reads, total 16S reads, total ARG reads, and total MGE reads
# This shell script recalls and requires "summary_reads_16S_ARG_MGE.py"
#------------------------------------------------------------------------------

####### Parameter setting ######################################################
# Your working directory 自分の作業ディレクトリ
DIR_WORKING="/home/ryohonda/GlobalAMR"

## the directory of raw sequence files (specify the absolute path. do not include the final '/')  
## 配列データのあるディレクトリ（絶対パスで指定。最後のスラッシュ '/' は含めない）
DIR_SEQ="/home/ryohonda/sequence"
DIR_QT="${DIR_SEQ}/1.qt" # directory of fastp reports
DIR_ARG="${DIR_WORKING}/4.ARG_profile" # directory of ARG profile of each sample
DIR_16S="${DIR_WORKING}/5.16S_taxonomy" # directory of 16S taxonomy of each sample
DIR_MGE="${DIR_WORKING}/7.MGE_profile" # directory of MGE profile of each sample

## output directory 
## 結果出力するディレクトリ（絶対パスで指定。最後のスラッシュ '/' は含めない）
DIR_SUM="${DIR_WORKING}/8.summary" # summary of reads

## location of the python script
PY_SUM="${DIR_WORKING}/scripts/summary_reads_16S_ARG_MGE.py"

####### Run program #############################################################
DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
echo "[$DATE] $JOB_NAME started. " 

## summarize reads data
python3 ${PY_SUM} ${DIR_QT} ${DIR_16S} ${DIR_ARG} ${DIR_MGE} ${DIR_SUM}

## record the completion of the script
DATE=$(date '+%Y-%m-%d %H:%M:%S %z')
echo "[$DATE] $JOB_NAME completed. " 
exit 0
