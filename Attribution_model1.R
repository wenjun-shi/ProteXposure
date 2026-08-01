########################
## Cox Attribution Analysis
########################

library(survival)
library(CoxR2)      
library(car)        
library(data.table)
library(optparse)

args_list = list(
  make_option(c("-d", "--disease"), type = "character", default="PHD09001",
              help="INPUT: disease code", metavar = "character"),
  make_option(c("-g", "--gender"), type = "character", default="female",
              help="INPUT: gender (male or female)", metavar = "character")
)
opt = parse_args(OptionParser(option_list=args_list))


PROJ_PATH <- "/public/home/gw_hychu/swj/"
DATA_PATH <- paste0(PROJ_PATH, "proteohubProject/02_data/")
OUT_PATH  <- paste0(PROJ_PATH, "web_out/Attribution_Cox1_", opt$gender, "/", opt$disease, "/") 

if (!dir.exists(OUT_PATH)) dir.create(OUT_PATH, recursive = TRUE)
suffix <- paste0("_", opt$gender)

disease_file <- paste0(DATA_PATH, "out_pheno_Cox/disease_trait", suffix, "_dat.rds")
follow_file  <- paste0(DATA_PATH, "follow_data/follow_data", suffix, ".rds")
cov_file     <- paste0(DATA_PATH, "cov/cov", suffix, ".rds")

disease_data <- readRDS(disease_file)
follow_data  <- readRDS(follow_file)
cov_data     <- readRDS(cov_file)

status <- as.numeric(disease_data[, opt$disease])
time   <- as.numeric(follow_data[[opt$disease]])

age_col <- if("Age" %in% colnames(cov_data)) "Age" else "age"
analysis_df <- data.frame(
  time = time,
  status = status,
  age = cov_data[[age_col]]
)

analysis_df <- na.omit(analysis_df)
analysis_df <- analysis_df[analysis_df$time > 0, ]

n_total <- nrow(analysis_df)
n_events <- sum(analysis_df$status)

## fit
surv_obj <- Surv(time = analysis_df$time, event = analysis_df$status)
model <- coxph(surv_obj ~ age, data = analysis_df)

res_r2_list <- coxr2(model)
final_cox_r2 <- as.numeric(res_r2_list$rsq)
c_index_val <- as.numeric(summary(model)$concordance)

# Wald Chi-square
anova_res <- car::Anova(model, type = "II", test = "Wald")
chisq_vals <- as.numeric(anova_res$`Chisq`)
chisq_percent <- (chisq_vals / sum(chisq_vals)) * 100
var_names <- as.character(rownames(anova_res))

results <- data.frame(
  Variable = var_names,
  Chisq = chisq_vals,
  Chisq_Contribution_Pct = chisq_percent,
  P_value = as.numeric(anova_res$`Pr(>Chisq)`),
  Model_CoxR2 = final_cox_r2,
  Model_C_Index = c_index_val,
  N = n_total,
  Events = n_events,
  stringsAsFactors = FALSE,
  row.names = NULL
)

output_file <- paste0(OUT_PATH, "cox1_", opt$gender, "_", opt$disease, ".tsv")
fwrite(results, file = output_file, sep = "\t")
