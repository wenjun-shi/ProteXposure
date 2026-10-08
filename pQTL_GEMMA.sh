#!/bin/bash
#SBATCH --job-name=mk_summ
#SBATCH --error=/public/home/gw_hychu/swj/geno/err_file/mk_summ1_%a.err
#SBATCH --out=/public/home/gw_hychu/swj/geno/out_file/mk_summ1_%a.out
#SBATCH --mem=2000
#SBATCH --array=1-3000%100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00

# Set parameters
GENO_PATH=/public/home/gw_hychu/swj/geno/EUR_protein/hm3/
PRO=/public/home/gw_hychu/swj/geno/04_olink/protein.txt
GEMMA=/public/home/gw_hychu/software/gemma-0.98.1-linux-static

k=${SLURM_ARRAY_TASK_ID}
pro=$(sed -n "${k}p" ${PRO})
if [ -z "${pro}" ]; then
    echo "No protein for task ID ${k}" >&2
    exit 1
fi

WORK_DIR="/public/home/gw_hychu/swj/geno/web_out/pqtl/${pro}"
for sex in male female; do
    INFILE="${GENO_PATH}${sex}/merge"
    ${GEMMA} -bfile ${INFILE} -n ${k} -notsnp -lm 1 -o summ_${sex}
    rm -f ./output/summ_${sex}.log.txt
    gzip -f ./output/summ_${sex}.assoc.txt
done