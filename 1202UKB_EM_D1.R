
#!/bin/bash

########################
## Estimate total effect for E2D and M2D
## Model_all: disease ~ exposure/measurement  + age + sex + BMI + site + 18 PCs
## Model_sex: disease ~ exposure/measurement  + age + BMI + site + 18 PCs

# Load packages
library(plyr)
library(dplyr)
library(data.table)
library(optparse)
library(ggplot2)

# Input parameters
args_list = list(
  make_option("--sex", type = "character", default=NULL,
              help="INPUT: sex label", metavar = "character"), 
  make_option("--fac_id", type = "character", default=NULL,
              help="INPUT: exposure factor id", metavar = "character")
)
opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

# opt = list(sex = "male", fac_id = "PHM01006")

# Set parameter
ALPHA <- 0.05
PROJ_PATH <- "/public/home/gw_hychu/swj"
DATA_PATH <- paste0(PROJ_PATH, "/proteohubProject/02_data")
dis_field_id <- fread(paste0(DATA_PATH, "/fliedID_disease_trait.txt"))

# Load data
dis_dat <- readRDS(paste0(DATA_PATH, "/out_pheno/disease_trait_", opt$sex, "_dat.rds"))

if(grepl("PHE", opt$fac_id)) {
  
  FAC_FILE <- paste0(DATA_PATH, "/out_pheno/exposure_factors_", opt$sex, ".rds")
  OUT_PATH <- paste0(PROJ_PATH, "/web_out/E-D/", opt$fac_id, "/")
  field_id <- fread(paste0(DATA_PATH, "/fieldID_exposure_factor.txt"))
}  
if(grepl("PHM", opt$fac_id)) {
  
  FAC_FILE <- paste0(DATA_PATH, "/out_pheno/measurement_trait_", opt$sex, ".rds")
  OUT_PATH <- paste0(PROJ_PATH, "/web_out/M-D/", opt$fac_id, "/")
  field_id <- fread(paste0(DATA_PATH, "/fieldID_measurement_trait.txt"))
}
fac_dat <- readRDS(FAC_FILE)
fac_dat_use <- fac_dat[, colnames(fac_dat) == opt$fac_id]
id_dat_all <- fread(paste0(DATA_PATH, "/sample_id/all_eid"))$V1
cov_dat <- readRDS(paste0(DATA_PATH, "/cov/cov_", opt$sex, ".rds"))
cov_dat$eid <- NULL
cov_dat$Array <- NULL
dummy <- model.matrix(~Site, data = cov_dat)
cov_dat <- cbind.data.frame(cov_dat, dummy[,-1])
cov_dat$Site <- NULL


# Fit total effect model
tot_summ <- alply(dis_dat, 2, function(pp){
  
  # pp=dis_dat[, 4]
  df_tmp <- data.frame(dis = pp, fac = fac_dat_use, cov_dat)
  
  if (all(is.na(pp))){
    
    cat("Sex specific trait.\n")
    res_tmp <- rep(NA, 4)
  } else{
    
    fit_mod <- glm(dis ~ ., data = df_tmp, family = "binomial")
    ci <- try(confint(fit_mod, names(coef(fit_mod))[2], method = "Wald") %>% round(5), silent = T)
    if (inherits(ci, "try-error")){
      
      res_tmp <- rep(NA, 4)
    } else {
      
      or_ci <- paste0(ci[1], ", ", ci[2])
      res_tmp <- c(coef(fit_mod)[2] , or_ci, coef(summary(fit_mod))[2, c(3, 4)])
    }
  }
  return(res_tmp)
  
}) %>% do.call("rbind", .)
total_result <- data.frame(FieldID = dis_field_id[[1]],
                           Trait = dis_field_id[[4]],
                           Category = dis_field_id[[2]],
                           beta = as.numeric(tot_summ[, 1]) %>% round(5),
                           CI =  tot_summ[, 2],
                           Z = tot_summ[, 3],
                           P = as.numeric(tot_summ[, 4]),
                           FactorID = opt$fac_id)
total_result <- na.omit(total_result)

# Figure
total_result_top <- total_result %>%
  mutate(low = as.numeric(tstrsplit(CI, ", ")[[1]])) %>%
  mutate(up = as.numeric(tstrsplit(CI, ", ")[[2]])) %>%
  arrange(P) %>%
  head(n=5)
total_result_top$Trait <- factor(total_result_top$Trait, levels = total_result_top$Trait)

# Plot 
plt <- ggplot(total_result_top, aes(x = Trait, y = beta, fill = Trait)) + 
  geom_bar(stat = "identity", width = 0.8) +
  geom_errorbar(aes(ymin = low, ymax = up), 
                size = 0.8, width = 0.2) +
  geom_text(aes(label = round(beta, 4)), vjust = -0.5,
            size = 4) +
  scale_fill_manual(values = c("#800000", "#2A9D8F", "#46F0F0", "#AAFFC3", "#6A3569")) +
  xlab("Top Five Diseases") +
  ylab("Effect Size") +
  ggtitle(paste0("Disease ~ ", field_id[[3]][field_id[[1]] == opt$fac_id], " + Covs"))+
  theme_bw() +
  theme(axis.text.x = element_text(size = 12, color = "black", angle = 70, hjust = 1),
        axis.title.x = element_text(size = 15, face = "bold", color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 15, face = "bold", color = "black"),
        panel.grid = element_blank(), 
        strip.placement = "outside",
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")

# Output
dir.create(OUT_PATH, recursive = TRUE, showWarnings = FALSE)
total_result <- total_result[total_result$P < ALPHA, ]
fwrite(total_result, sep = "\t",
       file = paste0(OUT_PATH, "/effect_size_", opt$sex, "_", opt$fac_id, ".tsv"))

ggsave(paste0(OUT_PATH, "/top_trait_", opt$sex, "_", opt$fac_id, ".png"),
       plt, 
       width = 6, height = 6, units = "in", dpi = 320)