########################
## Cox Attribution Analysis 3
########################

library(survival)
library(CoxR2)    
library(car)      
library(data.table)
library(optparse)

# 1. Input parameters
args_list = list(
  make_option("--disease", type = "character", default="PHD01001",
              help="INPUT: disease code", metavar = "character"),
  make_option("--imp", type = "integer", default = NULL,
              help="INPUT: imputation index (1-5)", metavar = "integer"),
  make_option("--gender", type = "character", default = "male",
              help="INPUT: gender (male or female)", metavar = "character")
)
opt = parse_args(OptionParser(option_list=args_list))

# Path configuration
PROJ_PATH <- "/public/home/gw_hychu/swj/"
DATA_PATH <- paste0(PROJ_PATH, "proteohubProject/02_data/")
SIG_PATH  <- paste0(PROJ_PATH, "web_out_Cox_", opt$gender, "/M-D/sig_PHM_", opt$gender, "/")
OUT_BASE  <- paste0(PROJ_PATH, "web_out/Attribution_Cox3_Gender/")

suffix_gender <- paste0("_", opt$gender)

## 2. Load Data
cat("--- Starting Analysis for", opt$gender, "| Disease:", opt$disease, "---\n")
disease_data <- readRDS(paste0(DATA_PATH, "out_pheno_Cox/disease_trait", suffix_gender, "_dat.rds"))
follow_data  <- readRDS(paste0(DATA_PATH, "follow_data/follow_data", suffix_gender, ".rds"))
cov_data     <- readRDS(paste0(DATA_PATH, "cov/cov", suffix_gender, ".rds"))
prs_all_dat  <- readRDS(paste0(DATA_PATH, "PRS", suffix_gender, ".rds"))

if (!is.null(opt$imp)) {
  cat("Using Imputed Dataset Index:", opt$imp, "\n")
  imp_file <- paste0(DATA_PATH, "out_pheno_attra/", opt$gender, "_imputed_PHM_list_5.rds")
  imputed_list <- readRDS(imp_file)
  phm_data <- as.matrix(imputed_list[[opt$imp]])
  rm(imputed_list); gc()
} else {
  cat("Warning: No imputation index specified. Using original raw data.\n")
  phm_data <- readRDS(paste0(DATA_PATH, "out_pheno_attra/measurement_trait_all.rds"))
}

sig_file <- paste0(SIG_PATH, opt$disease, "_", opt$gender, ".txt")
if (!file.exists(sig_file)) {
  stop(paste("No significant exposure file found at:", sig_file))
}

raw_content <- readLines(sig_file, warn = FALSE)
clean_vars <- unlist(strsplit(raw_content, "[, \t\n\r]+"))
clean_vars <- gsub('^c\\(|\\)$|^"|"$|^\'|\'$|,', '', clean_vars) 
sig_exposure_list <- clean_vars[clean_vars != "" & clean_vars != "c()"]


valid_exposures <- intersect(sig_exposure_list, colnames(phm_data))
status <- as.numeric(disease_data[, opt$disease])
time   <- as.numeric(follow_data[[opt$disease]])

prs_col_name <- paste0(opt$disease, "_Q")
has_prs <- prs_col_name %in% colnames(prs_all_dat)

if (has_prs) {
  cat("Having PRS \n")
  prs_vec <- as.numeric(as.character(prs_all_dat[[prs_col_name]]))
  pc_cols <- paste0("p22009_a", 1:18)
  
  analysis_df <- data.frame(
    time   = time,
    status = status,
    Age    = cov_data$Age,
    Array  = cov_data$Array,
    PRS_Q  = prs_vec,
    cov_data[, pc_cols, drop = FALSE]
  )

  covariate_names <- c(valid_exposures, "Age", "Array", "PRS_Q", pc_cols)
} else {
  analysis_df <- data.frame(
    time   = time,
    status = status,
    Age    = cov_data$Age
  )
  covariate_names <- c(valid_exposures, "Age")
}

exposure_matrix <- phm_data[, valid_exposures, drop = FALSE]
analysis_df <- cbind(analysis_df, as.data.frame(exposure_matrix))

analysis_df <- na.omit(analysis_df)
analysis_df <- analysis_df[analysis_df$time > 0, ]
n_total <- nrow(analysis_df)
n_events <- sum(analysis_df$status)

if(n_total < 30 | n_events < 5) stop("Insufficient samples or events.")

formula_str <- paste("Surv(time, status) ~", paste(covariate_names, collapse = " + "))
model <- coxph(as.formula(formula_str), data = analysis_df)

final_cox_r2 <- as.numeric(coxr2(model)$rsq)
c_index_val  <- as.numeric(summary(model)$concordance)

# Wald Chi-square
anova_res     <- car::Anova(model, type = "II", test = "Wald")
chisq_vals    <- as.numeric(anova_res$`Chisq`)
p_vals        <- as.numeric(anova_res$`Pr(>Chisq)`)
var_names     <- as.character(rownames(anova_res))

chisq_percent <- (chisq_vals / sum(chisq_vals)) * 100

results <- data.frame(
  Variable               = var_names,
  Chisq                  = chisq_vals,
  Chisq_Contribution_Pct = chisq_percent,
  P_value                = p_vals,
  stringsAsFactors       = FALSE
)

results$Model_CoxR2  <- final_cox_r2
results$Model_C_Index <- c_index_val
results$N            <- n_total
results$Events       <- n_events
results$Gender       <- opt$gender

OUT_PATH <- paste0(OUT_BASE, opt$gender, "/", opt$disease, "/")
if (!dir.exists(OUT_PATH)) dir.create(OUT_PATH, recursive = TRUE)

suffix_imp <- if(!is.null(opt$imp)) paste0("_imp", opt$imp) else ""
output_file <- paste0(OUT_PATH, "cox3_", opt$disease, suffix_imp, ".tsv")
fwrite(results, file = output_file, sep = "\t")

cat("\n=== Analysis Completed. Saved to:", output_file, "===\n")
