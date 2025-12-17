#!/bin/bash
#SBATCH --output=/home/username/log/%x_%j.out.log    # ★ユーザー名変更
#SBATCH --error=/home/username/log/%x_%j.err.log     # ★ユーザー名変更
#SBATCH -p medium
#SBATCH -N 1-1
#SBATCH -n 8
#SBATCH --mem=32G      # 100サンプル超えないなら十分(必要に応じて変更)

# スクリプト実行（★パス変更）
Rscript /home/username/16S_SILVA/scripts/s1_dada2_pipeline.r
