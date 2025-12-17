#!/bin/bash
#SBATCH --output=/home/ryohonda/log/%x_%j.out.log    # ★ユーザー名変更
#SBATCH --error=/home/ryohonda/log/%x_%j.err.log     # ★ユーザー名変更
#SBATCH -p short

Rscript /home/ryohonda/Sakai/s2_asv_postprocess.r
