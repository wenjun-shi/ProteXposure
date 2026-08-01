########################
## Cox Attribution Analysis 4
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
              help="INPUT: imputation index (1-5)", metavar = "integer")
)
opt = parse_args(OptionParser(option_list=args_list))

# Path configuration
PROJ_PATH <- "/public/home/gw_hychu/swj/"
DATA_PATH <- paste0(PROJ_PATH, "proteohubProject/02_data/")
SIG_PHM_PATH <- "/public/home/gw_hychu/swj/web_out_Cox_all/M-D/sig_PHM_all/"
SIG_PHE_PATH <- "/public/home/gw_hychu/swj/web_out_Cox_all/E-D/sig_PHE_all/"
OUT_BASE  <- paste0(PROJ_PATH, "web_out/Attribution_Cox4/")

## 2. Load Core Data
cat("Loading core datasets for:", opt$disease, "\n")
disease_data <- readRDS(paste0(DATA_PATH, "out_pheno_Cox/disease_trait_all_dat.rds"))
follow_data  <- readRDS(paste0(DATA_PATH, "follow_data/follow_data_all.rds"))
cov_data     <- readRDS(paste0(DATA_PATH, "cov/cov_all.rds"))
prs_all_dat  <- readRDS(paste0(DATA_PATH, "PRS_all.rds"))

if (!is.null(opt$imp)) {
  cat("Using Imputed Datasets (Index:", opt$imp, ")\n")

  phm_list <- readRDS(paste0(DATA_PATH, "out_pheno_attra/imputed_PHM_list_5.rds"))
  phm_data_mat <- as.matrix(phm_list[[opt$imp]])
  rm(phm_list); gc() 
  
  phe_list <- readRDS(paste0(DATA_PATH, "out_pheno_attra/imputed_PHE_list_5.rds"))
  phe_data_df <- phe_list[[opt$imp]]
  rm(phe_list); gc() 
  
} else {
  cat("Warning: No imputation index. Using original raw data.\n")
  phm_data_mat <- readRDS(paste0(DATA_PATH, "out_pheno_attra/measurement_trait_all.rds"))
  phe_data_df  <- readRDS(paste0(DATA_PATH, "out_pheno_attra/exposure_factors_all.rds"))
}

sig_phm_file <- paste0(SIG_PHM_PATH, opt$disease, "_all.txt")
valid_phm <- character(0)
if (file.exists(sig_phm_file)) {
  phm_vars <- as.character(fread(sig_phm_file)[[1]])
  valid_phm <- intersect(phm_vars, colnames(phm_data_mat))
}

sig_phe_file <- paste0(SIG_PHE_PATH, opt$disease, "_all.txt")
valid_phe <- character(0)
if (file.exists(sig_phe_file)) {
  phe_vars <- as.character(fread(sig_phe_file)[[1]])
  valid_phe <- intersect(phe_vars, colnames(phe_data_df))
}

all_valid_exposures <- c(valid_phm, valid_phe)
cat("Exposures to include: PHM(", length(valid_phm), "), PHE(", length(valid_phe), ")\n")

if(length(all_valid_exposures) == 0) {
  stop("No significant exposures found for this disease.")
}

status <- as.numeric(disease_data[, opt$disease])
time   <- as.numeric(follow_data[[opt$disease]])

prs_col_name <- paste0(opt$disease, "_Q")
has_prs <- prs_col_name %in% colnames(prs_all_dat)

if (has_prs) {
  prs_vec <- as.numeric(as.character(prs_all_dat[[prs_col_name]]))
  pc_cols <- paste0("p22009_a", 1:18)
  analysis_df <- data.frame(
    time   = time,
    status = status,
    Sex    = cov_data$Sex,
    Age    = cov_data$Age,
    Array  = cov_data$Array,
    PRS_Q  = prs_vec,
    cov_data[, pc_cols, drop = FALSE]
  )
  covariate_names <- c(all_valid_exposures, "Sex", "Age", "Array", "PRS_Q", pc_cols)
} else {
  analysis_df <- data.frame(
    time   = time,
    status = status,
    Sex    = cov_data$Sex,
    Age    = cov_data$Age
  )
  covariate_names <- c(all_valid_exposures, "Sex", "Age")
}

if(length(valid_phm) > 0) analysis_df <- cbind(analysis_df, as.data.frame(phm_data_mat[, valid_phm, drop = FALSE]))
if(length(valid_phe) > 0) analysis_df <- cbind(analysis_df, phe_data_df[, valid_phe, drop = FALSE])

analysis_df <- na.omit(analysis_df)
analysis_df <- analysis_df[analysis_df$time > 0, ]

n_total <- nrow(analysis_df); n_events <- sum(analysis_df$status)
if(n_total < 30 | n_events < 5) stop("Sample size or events too low.")

formula_str <- paste("Surv(time, status) ~", paste(covariate_names, collapse = " + "))
model <- coxph(as.formula(formula_str), data = analysis_df)

final_cox_r2 <- as.numeric(coxr2(model)$rsq)
c_index_val  <- as.numeric(summary(model)$concordance[1])
anova_res    <- car::Anova(model, type = "II", test = "Wald")

results <- data.frame(
  Variable                = rownames(anova_res),
  Chisq                  = as.numeric(anova_res$`Chisq`),
  Chisq_Contribution_Pct = (as.numeric(anova_res$`Chisq`) / sum(as.numeric(anova_res$`Chisq`))) * 100,
  P_value                = as.numeric(anova_res$`Pr(>Chisq)`),
  Model_CoxR2            = final_cox_r2,
  Model_C_Index          = c_index_val,
  N                      = n_total,
  Events                 = n_events,
  PRS_Included           = has_prs,
  stringsAsFactors       = FALSE
)


OUT_PATH <- paste0(OUT_BASE, opt$disease, "/")
if (!dir.exists(OUT_PATH)) dir.create(OUT_PATH, recursive = TRUE)

suffix <- if(!is.null(opt$imp)) paste0("_imp", opt$imp) else ""
fwrite(results, file = paste0(OUT_PATH, "cox4_", opt$disease, suffix, ".tsv"), sep = "\t")

cat("\n=== Analysis Completed. Saved with suffix:", suffix, "===\n")