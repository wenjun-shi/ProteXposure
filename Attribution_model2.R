########################
## Cox Attribution Analysis - Model 2 
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
OUT_PATH  <- paste0(PROJ_PATH, "web_out/Attribution_Cox2_", opt$gender, "/", opt$disease, "/") 

if (!dir.exists(OUT_PATH)) dir.create(OUT_PATH, recursive = TRUE)

disease_data <- readRDS(paste0(DATA_PATH, "out_pheno_Cox/disease_trait_", opt$gender, "_dat.rds"))
follow_data  <- readRDS(paste0(DATA_PATH, "follow_data/follow_data_", opt$gender, ".rds"))
cov_data     <- readRDS(paste0(DATA_PATH, "cov/cov_", opt$gender, ".rds"))
prs_all_dat  <- readRDS(paste0(DATA_PATH, "PRS_", opt$gender, ".rds"))

status <- as.numeric(disease_data[, opt$disease])
time   <- as.numeric(follow_data[, opt$disease])

# PRS
prs_col_name <- paste0(opt$disease, "_Q")
prs_vec <- as.numeric(as.character(prs_all_dat[[prs_col_name]]))
pc_cols <- paste0("p22009_a", 1:18)

analysis_df <- data.frame(
  time    = time,
  status  = status,
  Age     = as.numeric(cov_data[["Age"]]),
  Array   = as.factor(cov_data[["Array"]]),
  PRS_Q   = prs_vec
)

pc_matrix <- as.data.frame(lapply(cov_data[, pc_cols, drop = FALSE], as.numeric))
analysis_df <- cbind(analysis_df, pc_matrix)

analysis_df <- na.omit(analysis_df)
analysis_df <- analysis_df[analysis_df$time > 0, ]

n_total <- nrow(analysis_df)
n_events <- sum(analysis_df$status)

## fit
covariate_names <- c("Age", "Array", "PRS_Q", pc_cols)
formula_str <- paste("Surv(time, status) ~", paste(covariate_names, collapse = " + "))
model <- coxph(as.formula(formula_str), data = analysis_df)

res_r2_list <- coxr2(model)
final_cox_r2 <- as.numeric(res_r2_list$rsq)
c_index_val <- as.numeric(summary(model)$concordance)

anova_res <- car::Anova(model, type = "II", test = "Wald")
chisq_vals <- as.numeric(anova_res$`Chisq`)
chisq_percent <- (chisq_vals / sum(chisq_vals)) * 100
var_names <- as.character(rownames(anova_res))

n_vars <- length(var_names)
results <- data.frame(
  Variable               = var_names,
  Chisq                  = chisq_vals,
  Chisq_Contribution_Pct = chisq_percent,
  P_value                = as.numeric(anova_res$`Pr(>Chisq)`),
  Model_CoxR2            = rep(final_cox_r2, n_vars),
  Model_C_Index          = rep(c_index_val, n_vars),
  N                      = rep(n_total, n_vars),
  Events                 = rep(n_events, n_vars),
  stringsAsFactors       = FALSE
)

output_file <- paste0(OUT_PATH, "cox2_", opt$gender, "_", opt$disease, ".tsv")
fwrite(results, file = output_file, sep = "\t")
