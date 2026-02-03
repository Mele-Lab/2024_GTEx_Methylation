#!/bin/bash
#SBATCH --job-name=admixture
#SBATCH --chdir=/gpfs/scratch/bsc83/bsc83535/GTEx/v8/genotype_data
#SBATCH --output=/gpfs/scratch/bsc83/bsc83535/GTEx/v8/genotype_data/admixture2.out
#SBATCH --error=/gpfs/scratch/bsc83/bsc83535/GTEx/v8/genotype_data/admixture2.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=48 



for K in 3 4; 
do /gpfs/projects/bsc83/utils/admixture_linux-1.3.0/admixture --cv GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phasedLDPrunechr21.bed $K -j48  | tee log_LDPrunechr21_${K}.out;
done
