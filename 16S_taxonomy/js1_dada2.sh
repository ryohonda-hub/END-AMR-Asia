#!/bin/bash
#SBATCH --output=/home/ryohonda/log/%x_%j.out.log    # ★ユーザー名変更
#SBATCH --error=/home/ryohonda/log/%x_%j.err.log     # ★ユーザー名変更
#SBATCH -p medium
#SBATCH -N 1-1
#SBATCH -n 8
#SBATCH --mem=32G      # 45サンプルなら十分(必要に応じて変更)

# スクリプト実行（★パス変更）
Rscript /home/ryohonda/Sakai/s1_dada2_pipeline.r
