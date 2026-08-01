#!/bin/bash
#SBATCH -J PHE2D
#SBATCH -p batch
#SBATCH --array=41-54
#SBATCH --mem=20G
#SBATCH --time=1:00:00
#SBATCH -o /public/home/gw_hychu/BASE/PHE2D/std.out.%A_%a
#SBATCH -e /public/home/gw_hychu/BASE/PHE2D/std.err.%A_%a

source ~/.bashrc
conda activate my_r_env
cd /public/home/gw_hychu/swj

PHE_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" exposure_list.txt)
echo "task ${SLURM_ARRAY_TASK_ID}: ${PHE_ID}"

Rscript 1202UKB_EM_D.R \
  --sex all \
  --fac_id "${PHE_ID}" \
  > "logs/${PHE_ID}.log" 2>&1