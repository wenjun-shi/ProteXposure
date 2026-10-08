# Set parameters
GENO_PATH=/public/home/gw_hychu/swj/geno/EUR_protein/hm3/
DIS=/public/home/gw_hychu/swj/02_data/d_label.txt
GEMMA=/public/home/gw_hychu/software/gemma-0.98.1-linux-static
let k=0

for dis in $(cat ${DIS})
do
	let k=${k}+1
	if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]
	then
		for sex in male females
		# for sex in all
		do
			cd /public/home/gw_hychu/swj/web_out/GWAS-D/${dis}/
			INFILE=${GENO_PATH}${sex}/merge
			${GEMMA} -bfile ${INFILE} -n ${k} -notsnp -lm 1 -o summ_${sex}
			rm -rf ./output/summ_${sex}.log.txt
			gzip -f ./output/summ_${sex}.assoc.txt
		done
	fi
done
