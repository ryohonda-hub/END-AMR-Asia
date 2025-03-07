#!/bin/bash
#SBATCH --output=$HOME/log/%x_%j.out
#SBATCH --error=$HOME/log/%x_%j.err
#SBATCH -N 1-1 
#SBATCH -n 8

srun s3_16S-summary_ARG.sh
