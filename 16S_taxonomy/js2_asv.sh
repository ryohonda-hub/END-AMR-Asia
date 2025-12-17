#!/bin/bash
#SBATCH --output=/home/username/log/%x_%j.out.log    # ★ユーザー名変更
#SBATCH --error=/home/username/log/%x_%j.err.log     # ★ユーザー名変更
#SBATCH -p short

Rscript /home/username/16S_SILVA/scripts/s2_asv_postprocess.r
