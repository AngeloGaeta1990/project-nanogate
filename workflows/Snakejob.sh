#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N Nanogate
#$ -cwd
#$ -pe sharedmem 2
#$ -l h_rt=48:00:00
#$ -M s1899269@ed.ac.uk
#$ -m beas
#$ -l h_vmem=2G

module load anaconda/5.0.1
source activate Nanogate

snakemake --cores 2 -s pipelines/Snakefile

source deactivate
module unload anaconda
