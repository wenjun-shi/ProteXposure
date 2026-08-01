########################
## Interaction analysis with Cox Proportional Hazards Model
## Model: Surv(time, status) ~ snp * protein + 6 covariates
## h(t) = h0(t) × exp(β₁×SNP + β₂×Protein + β₃×Interaction + Σβᵢ×Covariate)
########################

# Load packages
library(plyr)
library(dplyr)
library(data.table)
library(optparse)
library(survival)

# Input parameters
args_list = list(
  make_option("--disease", type = "character", default="PHD01001",
              help="INPUT: disease code (e.g PHD01001)", metavar = "character"),
  make_option("--sex", type = "character", default="all",
              help="INPUT: sex (all, male, female)", metavar = "character")
)
opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

# Set parameter paths
PROJ_PATH <- "/public/home/gw_hychu/swj/"
DATA_PATH <- paste0(PROJ_PATH, "proteohubProject/02_data/")

## 1. Load disease data
disease_data <- readRDS(paste0(DATA_PATH, "out_pheno_Cox/disease_trait_", opt$sex, "_dat.rds"))
disease_status <- disease_data[, opt$disease]

## 2. Load follow-up data 
follow_data <- readRDS(paste0(DATA_PATH, "follow_data/follow_data_", opt$sex, ".rds"))
follow_time <- follow_data[[opt$disease]]

## 3. Load SNP list for this disease
snp_list_file <- paste0(DATA_PATH, "sig_SNP/sig_SNP_", opt$sex, "/", opt$disease, "_", opt$sex, ".txt")
snp_list <- fread(snp_list_file, header = FALSE)
snp_names <- snp_list[[1]]

## 4. Load genotype data
geno_mat <- fread(paste0(DATA_PATH, "sig_genotype/sig_geno_", opt$sex, ".raw"))
geno_colnames <- colnames(geno_mat)
geno_snp_ids <- laply(strsplit(geno_colnames, "_"), function(gn) gn[1])

# 提取可用的SNP数据
available_snps <- snp_names[snp_names %in% geno_snp_ids]
if(length(available_snps) == 0) {
  stop(paste("No available SNPs found for disease", opt$disease))
}

snp_data <- matrix(NA, nrow = nrow(geno_mat), ncol = length(available_snps))
colnames(snp_data) <- available_snps

for (i in 1:length(available_snps)) {
  snp_id <- available_snps[i]
  snp_idx <- which(geno_snp_ids == snp_id)
  if (length(snp_idx) > 1) {
    snp_idx <- snp_idx[1]
  }
  snp_data[, i] <- as.vector(geno_mat[[snp_idx]])
}

snp_data <- as.data.frame(snp_data)

## 5. Load significant proteins for this disease
sig_protein_file <- paste0(DATA_PATH, "sig_proteins_", opt$sex, "_Cox/", opt$disease, "_Proteins_", opt$sex, ".txt")
sig_proteins <- fread(sig_protein_file, header = FALSE)
pro_names <- sig_proteins[[1]]

## 6. Load protein data (only significant proteins)
all_protein_data <- readRDS(paste0(DATA_PATH, "04_olink/protein_clean_imp_", opt$sex, ".rds"))
pro_names <- intersect(pro_names, colnames(all_protein_data))
pro_data <- all_protein_data[, pro_names, drop = FALSE]

## 7. Load covariates
cov_data <- readRDS(paste0(DATA_PATH, "cov/cov_", opt$sex, ".rds"))
cov_data$eid <- NULL
dummy <- model.matrix(~Site, data = cov_data)
cov_complete <- cbind.data.frame(cov_data, dummy[,-1])
cov_complete$Site <- NULL

# Create output directory
BASE_OUT_PATH <- paste0(PROJ_PATH, "web_out_Cox/GP-D/single_results/")
disease_folder <- paste0(opt$disease, "_", opt$sex)
OUT_PATH <- paste0(BASE_OUT_PATH, disease_folder, "/")
if (!dir.exists(OUT_PATH)) {
  dir.create(OUT_PATH, recursive = TRUE)
}

# Create survival object
surv_obj <- Surv(time = follow_time, event = disease_status)

# Extract coefficients from Cox model and calculate HR
extract_cox_coefficients <- function(model, term_name) {
  coef_summary <- summary(model)$coefficients
  
  if(term_name %in% rownames(coef_summary)) {
    coef_row <- coef_summary[term_name, ]
    beta <- as.numeric(coef_row["coef"])
    se <- as.numeric(coef_row["se(coef)"])
    p_value <- as.numeric(coef_row["Pr(>|z|)"])
    
    z_value <- qnorm(0.975)
    hr <- exp(beta)
    hr_lower <- exp(beta - z_value * se)
    hr_upper <- exp(beta + z_value * se)
    
    return(c(beta, se, hr, hr_lower, hr_upper, p_value))
  } else {
    return(rep(NA, 6))
  }
}

# Single SNP interaction analysis using Cox model
all_results <- alply(available_snps, 1, function(snp_var){
  
  # Create base dataframe with survival object, SNP and covariates
  base_df <- data.frame(
    surv_obj = surv_obj,
    snp = snp_data[[snp_var]],
    cov_complete
  )
  
  # Analyze interaction with each significant protein
  protein_results <- alply(1:length(pro_names), 1, function(pro_idx){ 
    if(pro_idx %% 100 == 0) cat("  Protein:", pro_idx, "/", length(pro_names), "\n")
    tryCatch({
      # Add protein data
      analysis_df <- base_df
      analysis_df$protein <- pro_data[, pro_names[pro_idx]]
      analysis_df <- na.omit(analysis_df)    # Delete missing values before analysis
      
      if(sum(analysis_df$surv_obj[, 2]) < 10) {   # Check if we have enough events
        return(NA)
      }
      
      # Fit Cox interaction model with COMPLETE covariates     
      model <- coxph(surv_obj ~ snp * protein + ., 
                     data = analysis_df)
      
      snp_coef <- extract_cox_coefficients(model, "snp")        # Extract coefficients
      protein_coef <- extract_cox_coefficients(model, "protein")
      interaction_coef <- extract_cox_coefficients(model, "snp:protein")
      
      if(all(is.na(interaction_coef))) {
        interaction_coef <- extract_cox_coefficients(model, "protein:snp")
      }
      
      # result_vector
      result_vector <- c(
        as.character(snp_var),
        as.character(pro_names[pro_idx]),
        snp_coef,
        protein_coef,
        interaction_coef,
        nrow(analysis_df),
        sum(analysis_df$surv_obj[, 2])
      )
      
      return(result_vector)
      
    }, error = function(e) {
      return(NA)
    })
    
  })
  # Remove failed analyses
  protein_results <- protein_results[!is.na(protein_results)]
  if(length(protein_results) == 0) return(NULL)
  
  # Combine results for this SNP
  result_matrix <- do.call("rbind", protein_results)
  return(result_matrix)
})

# Remove NULL results
all_results <- all_results[!sapply(all_results, is.null)]

if(length(all_results) > 0) {
  final_results <- as.data.frame(do.call("rbind", all_results), stringsAsFactors = FALSE)
  
  colnames(final_results) <- c("SNP", "Protein", 
                               "Beta_SNP", "SE_SNP", "HR_SNP", 
                               "HR_Lower_95CI_SNP", "HR_Upper_95CI_SNP", "P_SNP",
                               "Beta_Protein", "SE_Protein", "HR_Protein", 
                               "HR_Lower_95CI_Protein", "HR_Upper_95CI_Protein", "P_Protein",
                               "Beta_Interaction", "SE_Interaction", "HR_Interaction", 
                               "HR_Lower_95CI_Interaction", "HR_Upper_95CI_Interaction", "P_Interaction",
                               "N", "Events")
  
  # Convert to the data type
  numeric_cols <- c("Beta_SNP", "SE_SNP", "HR_SNP", 
                    "HR_Lower_95CI_SNP", "HR_Upper_95CI_SNP", "P_SNP",
                    "Beta_Protein", "SE_Protein", "HR_Protein", 
                    "HR_Lower_95CI_Protein", "HR_Upper_95CI_Protein", "P_Protein",
                    "Beta_Interaction", "SE_Interaction", "HR_Interaction", 
                    "HR_Lower_95CI_Interaction", "HR_Upper_95CI_Interaction", "P_Interaction",
                    "N", "Events")
  
  for(col in numeric_cols) {
    final_results[[col]] <- as.numeric(final_results[[col]])
  }
  # Multiple testing correction (optional, added as in your logit code)
  final_results$P_Interaction_FDR <- p.adjust(final_results$P_Interaction, method = "fdr")
  final_results$P_Interaction_Bonferroni <- p.adjust(final_results$P_Interaction, method = "bonferroni")
  
  n_significant_05 <- sum(final_results$P_Interaction < 0.05, na.rm = TRUE)
  n_significant_fdr_05 <- sum(final_results$P_Interaction_FDR < 0.05, na.rm = TRUE)
  
  # Save all results
  output_file_full <- paste0(OUT_PATH, "cox_", opt$disease, "_", opt$sex, ".tsv")
  fwrite(final_results, file = output_file_full, sep = "\t")
  
  # Save significant results
  sig_results <- final_results[final_results$P_Interaction < 0.05, ]
  if(nrow(sig_results) > 0) {
    output_file_sig <- paste0(OUT_PATH, "cox_", opt$disease, "_", opt$sex, "_significant.tsv")
    fwrite(sig_results, file = output_file_sig, sep = "\t")
  }
  
  warnings()
  # Summary
  cat("\n=== Cox Model Analysis Summary ===\n")
  cat("Disease analyzed:", opt$disease, "\n")
  cat("Available SNPs for analysis:", length(available_snps), "\n")
  cat("Total interactions tested:", nrow(final_results), "\n")
  cat("Significant interactions (P_Interaction < 0.05):", n_significant_05, "\n")
  cat("FDR significant interactions (FDR < 0.05):", n_significant_fdr_05, "\n")
  
} else {
  cat("No successful analyses completed for disease", opt$disease, "\n")
}
