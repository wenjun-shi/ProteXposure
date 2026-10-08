library(bigreadr)
library(dplyr)
library(plyr)

# Set parameters
PROJ_PATH <- "/public/home/gw_hychu/swj/03_trait/proteohubProject"
TRAIT <- "/public/home/gw_hychu/swj/03_trait/participant.csv.gz"

# Load sample id and field id
sample_id <- list.files(paste0(PROJ_PATH, "02_data/sample_id/")) %>% 
  alply(., 1, function(ff) fread2(paste0(PROJ_PATH, "02_data/sample_id/", ff))[, 1])
id_name <- c("all", "female", "male")
field_id <- fread2(paste0(PROJ_PATH, "/02_data/fieldID_measurement_trait.txt"))

# Load traits
trait_eid <- fread2(TRAIT, select = "eid")
trait_quant_dat <- field_id$Code[!grepl("\\/", field_id$Code)] %>%
  fread2(TRAIT, select = .)
trait_quant_dat_ext <- field_id$Code[grepl("\\/", field_id$Code)] %>% 
  strsplit(., "\\/") %>%
  llply(., function(pp) {
    
    dat_tmp <- trait_quant_dat[, match(pp, colnames(trait_quant_dat))]
    return(dat_tmp[, 1] / dat_tmp[, 2])
  }) %>% do.call("cbind", .) 
colnames(trait_quant_dat_ext) <- field_id$Code[grepl("\\/", field_id$Code)] 
trait_quant_dat <- cbind.data.frame(trait_quant_dat, trait_quant_dat_ext)

# Output traits
trait_quant_sample <- alply(c(1: length(id_name)), 1, function (ss){
  
  trait_quant_dat_pro <- trait_quant_dat[match(sample_id[[ss]], trait_eid[, 1]), ] %>%
    apply(., 2, function(pp){ as.numeric(pp) }) 
  colnames(trait_quant_dat_pro) <- field_id$PHID[match(colnames(trait_quant_dat_pro), field_id$Code)]
  trait_quant_dat_pro <- trait_quant_dat_pro[, order(colnames(trait_quant_dat_pro))]
  trait_quant_summ <- apply(trait_quant_dat_pro, 2, function(tt) sum(!is.na(tt))) 
  saveRDS(trait_quant_dat_pro, file = paste0(PROJ_PATH, "02_data/out_pheno/measurement_trait_", id_name[ss], ".rds"))
  return(trait_quant_summ)
}) %>% do.call("cbind", .) 

write.csv(trait_quant_sample, quote = F,
          file = paste0(PROJ_PATH, "02_data/out_pheno/summary_measurement_trait.csv"))


