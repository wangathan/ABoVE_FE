#!/bin/bash
#SBATCH --time=9:00:00
#SBATCH --mem=100000
#SBATCH --mincpus=28

Rscript whichOnesArentDone.R $1 $2

