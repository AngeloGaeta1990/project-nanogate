#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N fast_test
#$ -cwd
#$ -pe smp 5
#$ -l h_rt=120:00:00
#$ -M s1899269@ed.ac.uk
#$ -m beas
#$ -l h_vmem=32G

source /export/homes/home/s1899269/miniconda3/bin/activate project-nanogate  
echo "source activated"

snakemake --cores 5 --configfile fast_test.yaml

conda deactivate


