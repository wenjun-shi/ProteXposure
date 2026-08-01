#!/bin/bash
#SBATCH -J PHM2D
#SBATCH -p batch
#SBATCH --array=123-129      
#SBATCH --mem=20G
#SBATCH --time=1:00:00
#SBATCH -o /public/home/gw_hychu/BASE/PHM2D/std.out.%A_%a
#SBATCH -e /public/home/gw_hychu/BASE/PHM2D/std.err.%A_%a

source ~/.bashrc
conda activate my_r_env
cd /public/home/gw_hychu/swj

PHM_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" measurement_list.txt)

echo "task ${SLURM_ARRAY_TASK_ID}: ${PHM_ID}"

Rscript 1202UKB_EM_D.R \
  --sex all \
  --fac_id "${PHM_ID}" \
  > "logs/${PHM_ID}.log" 2>&1