#!/bin/bash
#SBATCH --output=$HOME/log/%x_%j.out
#SBATCH --error=$HOME/log/%x_%j.err
#SBATCH -N 1-1 
#SBATCH -n 24
#SBATCH -t 9-23:59:59

srun s4_blastn-MGE_profile.sh
