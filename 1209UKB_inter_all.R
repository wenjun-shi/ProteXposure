########################
## Interaction analysis
## Model: disease ~ exposure * protein + 6 covariates
## Logit (p) = β₀+β₁× (M,E)+ β₂× (Protein)+β3× (Interaction)+ 6 Covariate
########################

library(plyr)
library(dplyr)
library(data.table)
library(optparse)

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
cat("Loading disease data for:", opt$disease, "\n")
t2d_data <- readRDS(paste0(DATA_PATH, "out_pheno/disease_trait_all_dat.rds"))

# Check if disease exists
if(!opt$disease %in% colnames(t2d_data)) {
  stop("Disease '", opt$disease, "' not found in data. Available diseases: ",
       paste(head(colnames(t2d_data), 5), collapse = ", "), "...")
}

disease_outcome <- t2d_data[, opt$disease]
cat("Successfully loaded disease outcome:", opt$disease, "\n")

## 2. Load exposure data based on prefix
cat("Loading exposure data for:", opt$exposure, "\n")

if(grepl("PHE", opt$exposure)) {
  # Exposure factor type
  exp_data <- readRDS(paste0(DATA_PATH, "out_pheno/exposure_factors_all.rds"))
  BASE_OUT_PATH <- paste0(PROJ_PATH, "web_out/EP-D/single_results/")
  cat("Exposure type: Exposure factor (PHE)\n")
} else if(grepl("PHM", opt$exposure)) {
  # Measurement trait type
  exp_data <- readRDS(paste0(DATA_PATH, "out_pheno/measurement_trait_all.rds"))
  BASE_OUT_PATH <- paste0(PROJ_PATH, "web_out/MP-D/single_results/")
  cat("Exposure type: Measurement trait (PHM)\n")
} else {
  stop("Exposure name '", opt$exposure, "' must start with PHE or PHM")
}

# OUT_PATH
disease_folder <- paste0(opt$disease, "_all")
OUT_PATH <- paste0(BASE_OUT_PATH, disease_folder, "/")

# Check if exposure exists
if(!opt$exposure %in% colnames(exp_data)) {
  stop("Exposure '", opt$exposure, "' not found in data. Available exposures: ",
       paste(head(colnames(exp_data), 5), collapse = ", "), "...")
}

## 3. Load significant proteins for this disease
cat("Loading significant proteins for disease:", opt$disease, "\n")
sig_protein_file <- paste0(DATA_PATH, "sig_proteins_all/", opt$disease, "_Proteins_all.txt")

if(!file.exists(sig_protein_file)) {
  stop("Significant protein file not found: ", sig_protein_file, 
       "\nPlease make sure the file exists in the correct format.")
}

# Read significant protein list
sig_proteins <- fread(sig_protein_file, header = FALSE)
if(ncol(sig_proteins) == 0 || nrow(sig_proteins) == 0) {
  stop("No significant proteins found in file: ", sig_protein_file)
}

# Get protein names - assuming first column is protein ID
pro_names <- sig_proteins[[1]]
cat("Found", length(pro_names), "significant proteins for disease", opt$disease, "\n")

## 4. Load protein data (only significant proteins)
all_protein_data <- readRDS(paste0(DATA_PATH, "04_olink/protein_clean_imp_all.rds"))

# Check if significant proteins exist in protein data
missing_proteins <- setdiff(pro_names, colnames(all_protein_data))
if(length(missing_proteins) > 0) {
  cat("Warning:", length(missing_proteins), "proteins not found in protein data:\n")
  print(head(missing_proteins))
  pro_names <- intersect(pro_names, colnames(all_protein_data))
  cat("Proceeding with", length(pro_names), "available proteins\n")
}

if(length(pro_names) == 0) {
  stop("No significant proteins available in the protein data for analysis.")
}

# Only select significant proteins
pro_data <- all_protein_data[, pro_names, drop = FALSE]
cat("Protein data dimensions:", dim(pro_data), "\n")

## 5. Load covariates
cov_data <- readRDS(paste0(DATA_PATH, "cov/cov_all.rds"))
cov_data$eid <- NULL
dummy <- model.matrix(~Site, data = cov_data)
cov_complete <- cbind.data.frame(cov_data, dummy[,-1])
cov_complete$Site <- NULL

cat("Using COMPLETE covariates:", colnames(cov_complete), "\n")
cat("Number of covariates:", ncol(cov_complete), "\n")

# Create output directory
cat("Creating output directory:", OUT_PATH, "\n")
if (!dir.exists(OUT_PATH)) {
  dir.create(OUT_PATH, recursive = TRUE)
}

# Select exposure to analyze
exp_names <- colnames(exp_data)
if(!is.null(opt$exposure)) {
  # Analyze a single specified exposure
  exp_batch <- opt$exposure
  if(!exp_batch %in% exp_names) {
    stop("Exposure '", exp_batch, "' not found in data. Available exposures: ", 
         paste(head(exp_names, 5), collapse = ", "), "...")
  }
  cat("Analyzing single exposure:", exp_batch, "with", length(pro_names), 
      "significant proteins for disease", opt$disease, "\n")
} else {
  # If not specified, analyze all exposures
  exp_batch <- exp_names
  cat("Analyzing all", length(exp_batch), "exposures with", 
      length(pro_names), "significant proteins for disease", opt$disease, "\n")
}

# Extract coefficients and calculate OR
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

# Single exposure interaction analysis
all_results <- alply(exp_batch, 1, function(exp_var){
  cat("Processing exposure:", exp_var, "\n")
  
  base_df <- data.frame(
    disease = disease_outcome,
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
      analysis_df <- na.omit(analysis_df)  # Delete missing values before analysis
      
      # Fit interaction model with COMPLETE covariates     
      model <- glm(disease ~ exposure * protein + ., 
                   data = analysis_df, 
                   family = binomial(),
                   control = glm.control(maxit = 50))
      # Extract coefficients
      exposure_coef <- extract_coefficients(model, "exposure")
      protein_coef <- extract_coefficients(model, "protein")
      interaction_coef <- extract_coefficients(model, "exposure:protein")
      
      if(all(is.na(interaction_coef))) {
        interaction_coef <- extract_coefficients(model, "protein:exposure")
      }
      
      # result_vector
      result_vector <- c(
        as.character(exp_var),
        as.character(pro_names[pro_idx]),
        exposure_coef,
        protein_coef,
        interaction_coef,
        nrow(analysis_df),
        as.numeric(model$converged),
        model$iter
      )
      
      return(result_vector)
      
    }, error = function(e) {
      cat("  Error with protein", pro_names[pro_idx], ":", e$message, "\n")
      return(NA)
    })
    
  }, .progress = "text")
  
  # Remove failed analyses
  protein_results <- protein_results[!is.na(protein_results)]
  if(length(protein_results) == 0) return(NULL)
  
  # Combine results for exposure
  result_matrix <- do.call("rbind", protein_results)
  return(result_matrix)
  
}, .progress = "text")

# Remove NULL results
all_results <- all_results[!sapply(all_results, is.null)]

if(length(all_results) > 0) {
  final_results <- as.data.frame(do.call("rbind", all_results), stringsAsFactors = FALSE)
  
  colnames(final_results) <- c("Exposure", "Protein", 
                               "Beta_Exposure", "SE_Exposure", "OR_Exposure", 
                               "OR_Lower_95CI_Exposure", "OR_Upper_95CI_Exposure", "P_Exposure",
                               "Beta_Protein", "SE_Protein", "OR_Protein", 
                               "OR_Lower_95CI_Protein", "OR_Upper_95CI_Protein", "P_Protein",
                               "Beta_Interaction", "SE_Interaction", "OR_Interaction", 
                               "OR_Lower_95CI_Interaction", "OR_Upper_95CI_Interaction", "P_Interaction",
                               "N", "Converged", "Iterations")
  
  # Convert to the data type
  numeric_cols <- c("Beta_Exposure", "SE_Exposure", "OR_Exposure", 
                    "OR_Lower_95CI_Exposure", "OR_Upper_95CI_Exposure", "P_Exposure",
                    "Beta_Protein", "SE_Protein", "OR_Protein", 
                    "OR_Lower_95CI_Protein", "OR_Upper_95CI_Protein", "P_Protein",
                    "Beta_Interaction", "SE_Interaction", "OR_Interaction", 
                    "OR_Lower_95CI_Interaction", "OR_Upper_95CI_Interaction", "P_Interaction",
                    "N", "Iterations")
  
  for(col in numeric_cols) {
    final_results[[col]] <- as.numeric(final_results[[col]])
  }
  
  final_results$Converged <- as.logical(as.numeric(final_results$Converged))
  
  # cheak data
  cat("Data preview:\n")
  print(head(final_results[, 1:8], 3))
  cat("\nColumn types:\n")
  print(sapply(final_results, class))
  
  n_significant <- sum(final_results$P_Interaction < 0.05, na.rm = TRUE)
  
  # Save all results
  if(!is.null(opt$exposure)) {
    output_file <- paste0(OUT_PATH, "interaction_", opt$exposure, "_all.tsv")
  } else {
    output_file <- paste0(OUT_PATH, "interactions_all.tsv")
  }
  
  fwrite(final_results, file = output_file, sep = "\t")
  