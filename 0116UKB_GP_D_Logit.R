########################
## model: disease ~ snp * protein + 6 covariates
## Logit (p) = β₀ + β₁× (SNP) + β₂× (Protein) + β₃× (SNP×Protein) + Covariates
########################

# Load packages
library(plyr)
library(dplyr)
library(data.table)
library(optparse)

# Input parameters
args_list = list(
  make_option("--disease", type = "character", default="PHD01001",
              help="INPUT: disease code (e.g PHD01001)", metavar = "character"),
  make_option("--sex", type = "character", default="all",
              help="INPUT: sex group (all, male, female)", metavar = "character")
)
opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

sex <- opt$sex
# Set parameter paths
PROJ_PATH <- "/public/home/gw_hychu/swj/"
DATA_PATH <- paste0(PROJ_PATH, "proteohubProject/02_data/")

## 1. Load disease data
t2d_data <- readRDS(paste0(DATA_PATH, "out_pheno/disease_trait_", sex, "_dat.rds"))
disease_outcome <- t2d_data[, opt$disease]

## 2. Load SNP list for this disease
snp_list_file <- paste0(DATA_PATH, "sig_SNP/sig_SNP_", sex, "/", opt$disease, "_", sex, ".txt")
snp_list <- fread(snp_list_file, header = FALSE)
snp_names <- snp_list[[1]]

## 3. Load genotype data
geno_mat <- fread(paste0(DATA_PATH, "sig_genotype/sig_geno_", sex, ".raw"))
geno_colnames <- colnames(geno_mat)
geno_snp_ids <- laply(strsplit(geno_colnames, "_"), function(gn) gn[1])

available_snps <- snp_names[snp_names %in% geno_snp_ids]
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

## 4. Load significant proteins for this disease
sig_protein_file <- paste0(DATA_PATH, "sig_proteins_", sex, "/", opt$disease, "_Proteins_", sex, ".txt")
sig_proteins <- fread(sig_protein_file, header = FALSE)
pro_names <- sig_proteins[[1]]

## 5. Load protein data (only significant proteins)
all_protein_data <- readRDS(paste0(DATA_PATH, "04_olink/protein_clean_imp_", sex, ".rds"))
pro_data <- all_protein_data[, pro_names, drop = FALSE]

## 6. Load covariates
cov_data <- readRDS(paste0(DATA_PATH, "cov/cov_", sex, ".rds"))
cov_data$eid <- NULL
dummy <- model.matrix(~Site, data = cov_data)
cov_complete <- cbind.data.frame(cov_data, dummy[,-1])
cov_complete$Site <- NULL

## 7. Prepare output directory
BASE_OUT_PATH <- paste0(PROJ_PATH, "web_out/GP-D/single_results/")
disease_folder <- paste0(opt$disease, "_", sex)
OUT_PATH <- paste0(BASE_OUT_PATH, disease_folder, "/")

if (!dir.exists(OUT_PATH)) {
  dir.create(OUT_PATH, recursive = TRUE)
}

## 8. Coefficient extraction function
extract_coefficients <- function(model, term_name) {
  coef_summary <- coef(summary(model))
  
  if(term_name %in% rownames(coef_summary)) {
    coef_row <- coef_summary[term_name, ]
    beta <- as.numeric(coef_row["Estimate"])
    se <- as.numeric(coef_row["Std. Error"])
    p_value <- as.numeric(coef_row["Pr(>|z|)"])
    z_value <- qnorm(0.975)
    or <- exp(beta)
    or_lower <- exp(beta - z_value * se)
    or_upper <- exp(beta + z_value * se)
    
    return(c(beta, se, or, or_lower, or_upper, p_value))
  } else {
    return(rep(NA, 6))
  }
}

## 9. Main analysis function
analyze_snp_protein_interaction <- function(snp_var, snp_data_mat, pro_data_mat, 
                                            disease_outcome, covariates) {
  base_df <- data.frame(
    disease = disease_outcome,
    snp = snp_data_mat[[snp_var]],
    covariates
  )
  
  protein_results <- alply(1:length(pro_names), 1, function(pro_idx) {
    tryCatch({
      analysis_df <- base_df
      analysis_df$protein <- pro_data_mat[, pro_names[pro_idx]]
      analysis_df <- na.omit(analysis_df)      
      
      # Fit interaction model
      model <- glm(disease ~ snp * protein + ., 
                   data = analysis_df, 
                   family = binomial(),
                   control = glm.control(maxit = 50))
      
      # Extract coefficients
      snp_coef <- extract_coefficients(model, "snp")
      protein_coef <- extract_coefficients(model, "protein")
      interaction_coef <- extract_coefficients(model, "snp:protein")
      if(all(is.na(interaction_coef))) {
        interaction_coef <- extract_coefficients(model, "protein:snp")
      }
      
      # Build result vector
      result_vector <- c(
        as.character(snp_var),
        as.character(pro_names[pro_idx]),
        snp_coef,
        protein_coef,
        interaction_coef,
        nrow(analysis_df)
      )
      
      return(result_vector)
      
    }, error = function(e) {
      return(NA)
    })
    
  })
  
  protein_results <- protein_results[!is.na(protein_results)]
  if(length(protein_results) == 0) {
    return(NULL)
  }
  
  result_matrix <- do.call("rbind", protein_results)
  return(result_matrix)
}

## 10. Execute main analysis
all_results <- alply(available_snps, 1, function(snp_var) {
  analyze_snp_protein_interaction(snp_var, snp_data, pro_data, 
                                  disease_outcome, cov_complete)
})

all_results <- all_results[!sapply(all_results, is.null)]

## 11. Process results
if(length(all_results) > 0) {
  final_results <- as.data.frame(do.call("rbind", all_results), stringsAsFactors = FALSE)
  
  colnames(final_results) <- c("SNP", "Protein", 
                               "Beta_SNP", "SE_SNP", "OR_SNP", 
                               "OR_Lower_95CI_SNP", "OR_Upper_95CI_SNP", "P_SNP",
                               "Beta_Protein", "SE_Protein", "OR_Protein", 
                               "OR_Lower_95CI_Protein", "OR_Upper_95CI_Protein", "P_Protein",
                               "Beta_Interaction", "SE_Interaction", "OR_Interaction", 
                               "OR_Lower_95CI_Interaction", "OR_Upper_95CI_Interaction", "P_Interaction",
                               "N")
  
  numeric_cols <- c("Beta_SNP", "SE_SNP", "OR_SNP", 
                    "OR_Lower_95CI_SNP", "OR_Upper_95CI_SNP", "P_SNP",
                    "Beta_Protein", "SE_Protein", "OR_Protein", 
                    "OR_Lower_95CI_Protein", "OR_Upper_95CI_Protein", "P_Protein",
                    "Beta_Interaction", "SE_Interaction", "OR_Interaction", 
                    "OR_Lower_95CI_Interaction", "OR_Upper_95CI_Interaction", "P_Interaction",
                    "N")
  
  for(col in numeric_cols) {
    final_results[[col]] <- as.numeric(final_results[[col]])
  }
  
  # Multiple testing correction
  final_results$P_Interaction_FDR <- p.adjust(final_results$P_Interaction, method = "fdr")
  final_results$P_Interaction_Bonferroni <- p.adjust(final_results$P_Interaction, method = "bonferroni")
  
  n_significant_05 <- sum(final_results$P_Interaction < 0.05, na.rm = TRUE)
  n_significant_fdr_05 <- sum(final_results$P_Interaction_FDR < 0.05, na.rm = TRUE)
  
  output_file_full <- paste0(OUT_PATH, "interaction_", opt$disease, "_", sex, ".tsv")
  fwrite(final_results, file = output_file_full, sep = "\t")
  
  sig_results <- final_results[final_results$P_Interaction < 0.05, ]
  if(nrow(sig_results) > 0) {
    output_file_sig <- paste0(OUT_PATH, "interaction_", opt$disease, "_", sex, "_significant.tsv")
    fwrite(sig_results, file = output_file_sig, sep = "\t")
  }
  
  cat(paste("Analysis completed:", opt$disease, "for sex =", sex, "\n"))
  cat(paste("Total interactions tested:", nrow(final_results), "\n"))
  cat(paste("Significant interactions (P < 0.05):", n_significant_05, "\n"))
  cat(paste("FDR significant interactions (FDR < 0.05):", n_significant_fdr_05, "\n"))
  
} else {
  cat("No successful analyses completed\n")
}