#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N similarity_07
#$ -cwd
#$ -pe smp 20
#$ -l h_rt=120:00:00
#$ -M s1899269@ed.ac.uk
#$ -m beas
#$ -l h_vmem=32G

source /export/homes/home/s1899269/miniconda3/bin/activate nanogate_env   
echo "source activated"

snakemake --cores 5 -s ../Snakefile --configfile similarity_07.yaml

conda deactivate


