########################
## Interaction analysis with Cox Proportional Hazards Model
## Model: Surv(time, status) ~ exposure * protein + 6 covariates
## h(t) = h0(t) × exp(β₁×Exposure + β₂×Protein + β₃×Interaction + Σβᵢ×Covariate)
########################

library(plyr)
library(dplyr)
library(data.table)
library(optparse)
library(survival)

# Input parameters
args_list = list(
  make_option("--exposure", type = "character", default=NULL, 
              help="INPUT: exposure variable name (PHE or PHM format)", metavar = "character"),
  make_option("--disease", type = "character", default="PHD03002",
              help="INPUT: disease code (e.g., PHD03002, PHD01001)", metavar = "character")
)
opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)
sex <- "all"

# Set parameter
PROJ_PATH <- "/public/home/gw_hychu/swj/"
DATA_PATH <- paste0(PROJ_PATH, "proteohubProject/02_data/")

## 1. Load disease data
diseas_data <- readRDS(paste0(DATA_PATH, "out_pheno_Cox/disease_trait_all_dat.rds"))
disease_status <- diseas_data[, opt$disease]

## 2. Load follow-up data 
follow_data <- readRDS(paste0(DATA_PATH, "follow_data/follow_data_all.rds"))
follow_time <- follow_data[[opt$disease]]

## 3. Load exposure data based on prefix
if(grepl("PHE", opt$exposure)) {
  exp_data <- readRDS(paste0(DATA_PATH, "out_pheno/exposure_factors_all.rds"))
  BASE_OUT_PATH <- paste0(PROJ_PATH, "web_out/EP-D/single_results/")
} else if(grepl("PHM", opt$exposure)) {
  exp_data <- readRDS(paste0(DATA_PATH, "out_pheno/measurement_trait_all.rds"))
  BASE_OUT_PATH <- paste0(PROJ_PATH, "web_out/MP-D/single_results/")
}

## 4. Load significant proteins for this disease
sig_protein_file <- paste0(DATA_PATH, "sig_proteins_all_Cox/", opt$disease, "_Proteins_all.txt")
sig_proteins <- fread(sig_protein_file, header = FALSE)
pro_names <- sig_proteins[[1]]

## 5. Load protein data (only significant proteins)
all_protein_data <- readRDS(paste0(DATA_PATH, "04_olink/protein_clean_imp_all.rds"))
pro_names <- intersect(pro_names, colnames(all_protein_data))
pro_data <- all_protein_data[, pro_names, drop = FALSE]

## 6. Load covariates
cov_data <- readRDS(paste0(DATA_PATH, "cov/cov_all.rds"))
cov_data$eid <- NULL
dummy <- model.matrix(~Site, data = cov_data)
cov_complete <- cbind.data.frame(cov_data, dummy[,-1])
cov_complete$Site <- NULL

# Create output directory
disease_folder <- paste0(opt$disease, "_all")
OUT_PATH <- paste0(BASE_OUT_PATH, disease_folder, "/")
if (!dir.exists(OUT_PATH)) {
  dir.create(OUT_PATH, recursive = TRUE)
}

# Select exposure to analyze
exp_batch <- opt$exposure

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

# Single exposure interaction analysis using Cox model
all_results <- alply(exp_batch, 1, function(exp_var){
  
  # Create base dataframe with survival object, exposure and covariates
  base_df <- data.frame(
    surv_obj = surv_obj,
    exposure = exp_data[[exp_var]],
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
      model <- coxph(surv_obj ~ exposure * protein + ., 
                     data = analysis_df)
      exposure_coef <- extract_cox_coefficients(model, "exposure")        # Extract coefficients
      protein_coef <- extract_cox_coefficients(model, "protein")
      interaction_coef <- extract_cox_coefficients(model, "exposure:protein")
      
      if(all(is.na(interaction_coef))) {
        interaction_coef <- extract_cox_coefficients(model, "protein:exposure")
      }
      
      # result_vector
      result_vector <- c(
        as.character(exp_var),
        as.character(pro_names[pro_idx]),
        exposure_coef,
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
  
  # Combine results for exposure
  result_matrix <- do.call("rbind", protein_results)
  return(result_matrix)
  
})

# Remove NULL results
all_results <- all_results[!sapply(all_results, is.null)]

if(length(all_results) > 0) {
  final_results <- as.data.frame(do.call("rbind", all_results), stringsAsFactors = FALSE)
  
  colnames(final_results) <- c("Exposure", "Protein", 
                               "Beta_Exposure", "SE_Exposure", "HR_Exposure", 
                               "HR_Lower_95CI_Exposure", "HR_Upper_95CI_Exposure", "P_Exposure",
                               "Beta_Protein", "SE_Protein", "HR_Protein", 
                               "HR_Lower_95CI_Protein", "HR_Upper_95CI_Protein", "P_Protein",
                               "Beta_Interaction", "SE_Interaction", "HR_Interaction", 
                               "HR_Lower_95CI_Interaction", "HR_Upper_95CI_Interaction", "P_Interaction",
                               "N", "Events")
  
  # Convert to the data type
  numeric_cols <- c("Beta_Exposure", "SE_Exposure", "HR_Exposure", 
                    "HR_Lower_95CI_Exposure", "HR_Upper_95CI_Exposure", "P_Exposure",
                    "Beta_Protein", "SE_Protein", "HR_Protein", 
                    "HR_Lower_95CI_Protein", "HR_Upper_95CI_Protein", "P_Protein",
                    "Beta_Interaction", "SE_Interaction", "HR_Interaction", 
                    "HR_Lower_95CI_Interaction", "HR_Upper_95CI_Interaction", "P_Interaction",
                    "N", "Events")
  
  for(col in numeric_cols) {
    final_results[[col]] <- as.numeric(final_results[[col]])
  }
  
  # Save all results
  output_file <- paste0(OUT_PATH, "cox_", opt$exposure, "_all.tsv")
  fwrite(final_results, file = output_file, sep = "\t")
