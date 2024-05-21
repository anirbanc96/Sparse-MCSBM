#!/bin/bash
#$ -N OverRatio           # <- name your job
#$ -j y              # <- join output and error for easy reading
#$ -pe mpi 25     # <- requesting 8 cores
#$ -l m_mem_free=5G  # <- start with MUCH LESS RAM

time Rscript --no-save BodyComp.R
