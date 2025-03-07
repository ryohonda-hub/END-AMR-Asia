#!/bin/bash
#SBATCH --output=$HOME/log/%x_%j.out
#SBATCH --error=$HOME/log/%x_%j.err

srun s5_summary_ARGMGE.sh
