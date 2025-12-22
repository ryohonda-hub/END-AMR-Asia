#!/bin/bash
#SBATCH --output=/home/username/log/%x_%j.out.log    # ★ユーザー名変更 / change to your username
#SBATCH --error=/home/username/log/%x_%j.err.log     # ★ユーザー名変更 / change to your username
#SBATCH -p short

# スクリプト実行（★パス変更）/ Change to the path of your script.
Rscript /home/username/16S_SILVA/scripts/s2_asv_postprocess.r
