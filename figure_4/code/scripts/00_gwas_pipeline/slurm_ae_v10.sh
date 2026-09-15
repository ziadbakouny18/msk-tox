#!/bin/bash
#SBATCH --job-name=ae_gwas_v10
#SBATCH --partition=componc_cpu
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --time=2-00:00:00
#SBATCH --output=/data1/reznike/walserr_gwas/logs/ae_gwas_v10_task%a_%j.out
#SBATCH --error=/data1/reznike/walserr_gwas/logs/ae_gwas_v10_task%a_%j.err
echo "v10 task ${SLURM_ARRAY_TASK_ID} started on $(hostname) at $(date)"
cd /data1/reznike/walserr_gwas
source activate gwas_env 2>/dev/null || conda activate gwas_env 2>/dev/null || true
Rscript 02_ae_gwas_worker_v10.R
echo "v10 task ${SLURM_ARRAY_TASK_ID} completed at $(date)"
