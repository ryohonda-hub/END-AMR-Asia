#!/bin/bash
#SBATCH --output=/home/username/log/%x_%j.out.log    # ★ユーザー名変更 / change to your username
#SBATCH --error=/home/username/log/%x_%j.err.log     # ★ユーザー名変更 / change to your username
#SBATCH -p medium
#SBATCH -N 1-1
#SBATCH -n 8
#SBATCH --mem=32G      # 100サンプル超えないなら32GBで十分(必要に応じて変更) / 32GB is enough for <100 samples

# スクリプト実行（★パス変更）/ Change to the path of your script.
Rscript /home/username/16S_SILVA/scripts/s1_dada2_pipeline.r
