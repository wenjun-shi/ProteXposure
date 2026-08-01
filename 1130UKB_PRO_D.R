########################
## Estimate total effect for Protein to Disease (P2D)
## Model: disease ~ protein + age + sex + BMI + site + 18 PCs

# Load packages
library(plyr)
library(dplyr)
library(data.table)
library(optparse)

# Input parameters
args_list = list(
  make_option("--sex", type = "character", default=NULL,
              help="INPUT: sex label", metavar = "character"), 
  make_option("--protein_list", type = "character", default=NULL,
              help="INPUT: protein list file", metavar = "character"),
  make_option("--disease_id", type = "character", default=NULL,
              help="INPUT: disease id to analyze (use 'all' for all diseases)", metavar = "character")
)
opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

# Set parameter
ALPHA <- 0.05
PROJ_PATH <- "/public/home/gw_hychu/swj"
DATA_PATH <- paste0(PROJ_PATH, "/proteohubProject/02_data")
dis_field_id <- fread(paste0(DATA_PATH, "/fliedID_disease_trait.txt"))

# Load data
dis_dat <- readRDS(paste0(DATA_PATH, "/out_pheno/disease_trait_", opt$sex, "_dat.rds"))

# Load protein data 
PROTEIN_FILE <- paste0(DATA_PATH, "/04_olink/protein_clean_imp_", opt$sex, ".rds")
protein_dat <- readRDS(PROTEIN_FILE)

# Load protein list
if (!is.null(opt$protein_list)) {
  protein_ids <- fread(opt$protein_list, header = FALSE)$V1
} else {
  protein_ids <- colnames(protein_dat)
}

# Select diseases to analyze
if (!is.null(opt$disease_id) && opt$disease_id != "all") {
  disease_ids <- strsplit(opt$disease_id, ",")[[1]]
  disease_ids <- disease_ids[disease_ids %in% colnames(dis_dat)]
  if (length(disease_ids) == 0) {
    stop("No valid disease IDs found in the data")
  }
} else {
  disease_ids <- colnames(dis_dat)
}

cat("Analyzing", length(disease_ids), "diseases\n")

id_dat_all <- fread(paste0(DATA_PATH, "/sample_id/all_eid"))$V1
cov_dat <- readRDS(paste0(DATA_PATH, "/cov/cov_", opt$sex, ".rds"))
cov_dat$eid <- NULL
dummy <- model.matrix(~Site, data = cov_dat)
cov_dat <- cbind.data.frame(cov_dat, dummy[,-1])
cov_dat$Site <- NULL

# Create output directory
OUT_PATH <- paste0(PROJ_PATH, "/web_out/P2D/")
dir.create(OUT_PATH, recursive = TRUE, showWarnings = FALSE)

# Initialize list to store results for each disease
disease_results <- list()

# Loop through each disease
for (disease_id in disease_ids) {
  
  cat("Analyzing disease:", disease_id, "\n")
  
  # Extract disease data
  disease_vector <- dis_dat[, disease_id]
  
  # Get disease information
  dis_info <- dis_field_id[dis_field_id[[1]] == disease_id, ]
  if (nrow(dis_info) == 0) {
    trait_name <- disease_id
    category_name <- "Unknown"
  } else {
    trait_name <- ifelse(ncol(dis_info) >= 4, dis_info[[4]], disease_id)
    category_name <- ifelse(ncol(dis_info) >= 2, dis_info[[2]], "Unknown")
  }
  
  # Initialize results for current disease
  disease_result <- data.frame()
  proteins_analyzed <- 0
  
  # Loop through all proteins for current disease
  for (protein_id in protein_ids) {
    
    # Check if protein exists in the data
    if (!protein_id %in% colnames(protein_dat)) {
      next
    }
    
    # Extract protein data
    if (is.matrix(protein_dat)) {
      protein_vector <- protein_dat[, protein_id]
    } else {
      protein_vector <- protein_dat[[protein_id]]
    }
    
    # Create data frame
    df_tmp <- data.frame(
      dis = disease_vector,
      protein = protein_vector,
      cov_dat
    )
    
    # Remove rows with missing
    complete_cases <- complete.cases(df_tmp)
    n_samples <- sum(complete_cases)
    df_complete <- df_tmp[complete_cases, ]
    
    if (n_samples == 0) {
      next
    } else if (length(unique(df_complete$dis)) < 2) {
      next
    } else {
      # Fit model
      fit_mod <- try(glm(dis ~ ., data = df_complete, family = "binomial"), silent = TRUE)
      
      if (inherits(fit_mod, "try-error")) {
        next
      }
      
      # Calculate confidence interval for beta
      ci <- try(confint(fit_mod, "protein", level = 0.95), silent = TRUE)
      
      if (inherits(ci, "try-error")){
        next
      } else {
        # Extract beta and its CI
        beta <- coef(fit_mod)["protein"]
        beta_lower <- ci[1]
        beta_upper <- ci[2]
        beta_ci <- paste0(sprintf("%.5e", beta_lower), ", ", sprintf("%.5e", beta_upper))
        
        # Extract statistics
        coef_summary <- coef(summary(fit_mod))
        if ("protein" %in% rownames(coef_summary)) {
          z_value <- coef_summary["protein", "z value"]
          p_value <- coef_summary["protein", "Pr(>|z|)"]
        } else {
          next
        }
        
        # Create result row
        protein_result <- data.frame(
          FieldID = disease_id,
          Trait = trait_name,
          Category = category_name,
          Protein = protein_id,
          N = n_samples,
          beta = sprintf("%.5e", beta),
          beta_CI = beta_ci,
          Z = sprintf("%.5f", z_value),
          P = sprintf("%.5e", p_value),
          stringsAsFactors = FALSE
        )
        
        disease_result <- rbind(disease_result, protein_result)
        proteins_analyzed <- proteins_analyzed + 1
        
        # Progress indicator
        if (proteins_analyzed %% 500 == 0) {
          cat("  Processed", proteins_analyzed, "proteins\n")
        }
      }
    }
  }
  
  # Perform multiple testing correction for current disease
  if (nrow(disease_result) > 0) {
    # Convert P values to numeric for correction
    disease_result$P_numeric <- as.numeric(disease_result$P)
    disease_result$P_Bon <- p.adjust(disease_result$P_numeric, method = "bonferroni")
    disease_result$P_BH <- p.adjust(disease_result$P_numeric, method = "BH")
    
    # Format adjusted P values
    disease_result$P_Bon <- sprintf("%.5e", disease_result$P_Bon)
    disease_result$P_BH <- sprintf("%.5e", disease_result$P_BH)
    
    # Remove temporary numeric column
    disease_result$P_numeric <- NULL
    
    # Sort by P value
    disease_result <- disease_result[order(as.numeric(disease_result$P)), ]
    
    # Store in list
    disease_results[[disease_id]] <- disease_result
    
    # Output file for current disease 
    disease_output_file <- paste0(OUT_PATH, "effect_size_", disease_id, "_", opt$sex, ".tsv")
    
    fwrite(disease_result, sep = "\t", file = disease_output_file)
    
    cat("  ", proteins_analyzed, "proteins analyzed,", 
        sum(as.numeric(disease_result$P_BH) < ALPHA), "FDR-significant\n")
  } else {
    cat("  No proteins successfully analyzed\n")
  }
}

# Save summary results
if (length(disease_results) > 0) {
  all_results <- do.call(rbind, disease_results)
  summary_file <- paste0(OUT_PATH, "P2D_summary_", opt$sex, ".tsv")
  fwrite(all_results, sep = "\t", file = summary_file)
  cat("\nCombined results saved to:", summary_file, "\n")
  
  # Final summary
  total_associations <- nrow(all_results)
  total_fdr_sig <- sum(as.numeric(all_results$P_BH) < ALPHA)
  total_bon_sig <- sum(as.numeric(all_results$P_Bon) < ALPHA)
  
  cat("\n=== Final Summary ===\n")
  cat("FDR-significant (BH<0.05):", total_fdr_sig, "\n")
}

cat("\nAnalysis completed\n")