#!/bin/bash
# Usage: sbatch slurm-serial-job-script
# Author Feng Yan, ACBD
#     feng.yan@monash.edu
# NOTE: To activate a SLURM option, remove the whitespace between the '#' and 'SBATCH'
# To give your job a name, replace "MyJob" with an appropriate name
#SBATCH --job-name=SA-atac
# To set a project account for credit charging, 
#SBATCH --account=ls25
# Request CPU resource for a serial job
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
# Memory usage (MB)
#SBATCH --mem-per-cpu=4000
# Set your minimum acceptable walltime, format: day-hours:minutes:seconds
#SBATCH --time=4:00:00
# To receive an email when job completes or fails
#SBATCH --mail-user=feng.yan@monash.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
# SBATCH --array=1-6
# Use reserved node to run job when a node reservation is made for you already
# SBATCH --reservation=reservation_name
#SBATCH --partition=genomics
#SBATCH --qos=genomics
# Command to run a serial job


module load anaconda/2019.03-Python3.7-gcc5
source activate /home/fyan0011/mx82/fyan0011/conda_env/rnasik-1.5.4

nextflow run alexyfyf/atac_nf --reads 'fastq/*R{1,2}.fastq.gz' -resume --fasta /home/fyan0011/ls25/feng.yan_raw/refFiles/ucsc_mm10/mm10.fa -profile slurm --outdir result_atac

