########################
## Estimate total effect for E2D and M2D using Cox proportional hazards model
## Model: Surv(time, status) ~ exposure/measurement + age + sex + BMI + site + 18 PCs

# Load packages
library(plyr)
library(dplyr)
library(data.table)
library(optparse)
library(ggplot2)
library(survival)  # For Cox models

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
dis_dat <- readRDS(paste0(DATA_PATH, "/out_pheno_Cox/disease_trait_", opt$sex, "_dat.rds"))  ##change

# Load follow-up time data
followup_dat <- readRDS(paste0(DATA_PATH, "/follow_data/follow_data_", opt$sex, ".rds"))

if(grepl("PHE", opt$fac_id)) {
  
  FAC_FILE <- paste0(DATA_PATH, "/out_pheno/exposure_factors_", opt$sex, ".rds")
  OUT_PATH <- paste0(PROJ_PATH, "/web_out_Cox/E-D/", opt$fac_id, "/")     ##change
  field_id <- fread(paste0(DATA_PATH, "/fieldID_exposure_factor.txt"))
}  
if(grepl("PHM", opt$fac_id)) {
  
  FAC_FILE <- paste0(DATA_PATH, "/out_pheno/measurement_trait_", opt$sex, ".rds")
  OUT_PATH <- paste0(PROJ_PATH, "/web_out_Cox/M-D/", opt$fac_id, "/")    ##change
  field_id <- fread(paste0(DATA_PATH, "/fieldID_measurement_trait.txt"))
}


# 创建输出目录（如果不存在）
if (!dir.exists(OUT_PATH)) {
  dir.create(OUT_PATH, recursive = TRUE, showWarnings = FALSE)
  cat("Created output directory:", OUT_PATH, "\n")
}

fac_dat <- readRDS(FAC_FILE)
fac_dat_use <- fac_dat[, colnames(fac_dat) == opt$fac_id]
id_dat_all <- fread(paste0(DATA_PATH, "/sample_id/all_eid"))$V1
cov_dat <- readRDS(paste0(DATA_PATH, "/cov/cov_", opt$sex, ".rds"))
cov_dat$eid <- NULL
dummy <- model.matrix(~Site, data = cov_dat)
cov_dat <- cbind.data.frame(cov_dat, dummy[,-1])
cov_dat$Site <- NULL

# Fit Cox proportional hazards model
tot_summ <- alply(1:ncol(dis_dat), 1, function(i){
  
  # Get disease-specific status and time
  disease_status <- dis_dat[, i]
  disease_time <- followup_dat[, i]
  
  # Create survival data frame
  df_tmp <- data.frame(
    time = disease_time,
    status = disease_status,
    fac = fac_dat_use,
    cov_dat
  )
  
  # Remove rows with missing values
  df_tmp <- na.omit(df_tmp)
  
  if (nrow(df_tmp) == 0){
    
    cat("Sex specific trait or no data.\n")
    res_tmp <- rep(NA, 5)
  } else if (sum(df_tmp$status == 1, na.rm = TRUE) < 5) {
    
    cat("Insufficient events (<5) for Cox model.\n")
    res_tmp <- rep(NA, 5)
  } else{
    
    tryCatch({
      # Fit Cox proportional hazards model
      fit_mod <- coxph(Surv(time, status) ~ ., data = df_tmp)
      
      # Extract coefficient and confidence interval for the factor
      beta <- coef(fit_mod)[1]
      hr <- exp(beta)
      ci <- exp(confint(fit_mod, parm = 1, level = 0.95))
      hr_ci <- paste0(round(ci[1], 5), ", ", round(ci[2], 5))
      
      # Get Z and P values
      coef_summary <- summary(fit_mod)
      z <- coef_summary$coefficients[1, "z"]
      p <- coef_summary$coefficients[1, "Pr(>|z|)"]
      
      res_tmp <- c(beta, hr_ci, z, p, hr)
    }, error = function(e) {
      cat("Error in Cox model:", e$message, "\n")
      res_tmp <- rep(NA, 5)
    })
  }
  return(res_tmp)
  
}) %>% do.call("rbind", .)

total_result <- data.frame(FieldID = dis_field_id[[1]],
                           Trait = dis_field_id[[4]],
                           Category = dis_field_id[[2]],
                           beta = as.numeric(tot_summ[, 1]) %>% round(5),
                           HR = as.numeric(tot_summ[, 5]) %>% round(5),  # Hazard Ratio
                           CI =  tot_summ[, 2],  # Now contains HR confidence interval
                           Z = tot_summ[, 3],
                           P = as.numeric(tot_summ[, 4]),
                           FactorID = opt$fac_id)
total_result <- na.omit(total_result)

# Figure - modified to show HR instead of beta
total_result_top <- total_result %>%
  mutate(low = as.numeric(tstrsplit(CI, ", ")[[1]])) %>%
  mutate(up = as.numeric(tstrsplit(CI, ", ")[[2]])) %>%
  arrange(P) %>%
  head(n=5)
total_result_top$Trait <- factor(total_result_top$Trait, levels = total_result_top$Trait)

# Plot - modified to show HR instead of beta
plt <- ggplot(total_result_top, aes(x = Trait, y = HR, fill = Trait)) + 
  geom_bar(stat = "identity", width = 0.8) +
  geom_errorbar(aes(ymin = low, ymax = up), 
                size = 0.8, width = 0.2) +
  geom_text(aes(label = round(HR, 4)), vjust = -0.5,
            size = 4) +
  scale_fill_manual(values = c("#800000", "#2A9D8F", "#46F0F0", "#AAFFC3", "#6A3569")) +
  xlab("Top Five Diseases") +
  ylab("Hazard Ratio") +
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
plt <- plt + expand_limits(y = max(total_result_top$up) * 1.03)
# Output
total_result <- total_result[total_result$P < ALPHA, ]
fwrite(total_result, sep = "\t",
       file = paste0(OUT_PATH, "/effect_size_", opt$sex, "_", opt$fac_id, ".tsv"))

ggsave(paste0(OUT_PATH, "/top_trait_", opt$sex, "_", opt$fac_id, ".png"),
       plt, 
       width = 6, height = 6, units = "in", dpi = 320)