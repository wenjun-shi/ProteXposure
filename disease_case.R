library(bigreadr)
library(plyr)
library(dplyr)

# Set parameters
PROJ_PATH <- "/public/home/gw_hychu/swj/03_trait/proteohubProject"
TRAIT <- "/public/home/gw_hychu/swj/03_trait/participant.csv.gz"
NONCANCER_DATE <- "/public/home/gw_hychu/swj/03_trait/noncancer_date.rds"
NONCANCER_ID <- "/public/home/gw_hychu/swj/03_trait/fieldID_noncancer_date.txt"

# Load sample id
all_id <- fread2(paste0(PROJ_PATH, "/sample_id/all_eid"))[, 1]
male_id <- fread2(paste0(PROJ_PATH, "/sample_id/male_eid"))[, 1]
female_id <- fread2(paste0(PROJ_PATH, "/sample_id/female_eid"))[, 1]

# Load disease
## Set field id
cancer_self_report <- c(paste0("p20001_i0_a", c(0:5)), paste0("p20001_i1_a", c(0:5)),
                        paste0("p20001_i2_a", c(0:5)), paste0("p20001_i3_a", c(0:5)))
cancer_self_report_age <- c(paste0("p20007_i0_a", c(0:5)), paste0("p20007_i1_a", c(0:5)),
                            paste0("p20007_i2_a", c(0:5)), paste0("p20007_i3_a", c(0:5)))
cancer_icd10 <- paste0("p40006_i", c(0: 21))
cancer_icd10_age <- paste0("p40008_i", c(0: 21))
noncancer_self_report <- c(paste0("p20002_i0_a", c(0:33)), paste0("p20002_i1_a", c(0:33)),
                           paste0("p20002_i2_a", c(0:33)), paste0("p20002_i3_a", c(0:33)))
noncancer_self_report_age <- c(paste0("p20009_i0_a", c(0:33)), paste0("p20009_i1_a", c(0:33)),
                               paste0("p20009_i2_a", c(0:33)), paste0("p20009_i3_a", c(0:33)))

## Load binary data
trait_binary_dat <- fread2(TRAIT, 
                           select = c("eid", 
                                      "p21022",         # Age at recruitment
                                      "p34",            # Year of birth
                                      "p41202",         # Main diagnose
                                      "p41204",         # Second diagnose
                                      cancer_self_report, cancer_self_report_age,
                                      cancer_icd10, cancer_icd10_age, 
                                      noncancer_self_report, noncancer_self_report_age))
trait_binary_dat <- trait_binary_dat[match(all_id, trait_binary_dat$eid), ]
sample_size <- nrow(trait_binary_dat)
recruitment_age <- trait_binary_dat$p21022
noncancer_date_df <- readRDS(NONCANCER_DATE)
noncancer_date_df <- noncancer_date_df[match(all_id, noncancer_date_df$eid), ]
noncancer_id_date <- fread2(NONCANCER_ID, header = F)

# Process field ID data
field_id <- fread2(paste0(PROJ_PATH, "/02_data/fliedID_disease_trait.txt"))
field_id$ICD10 <- gsub(",", "\\|", field_id$ICD10)

# Identify cancer and extract ages
trait_binary_dat_cancer <- trait_binary_dat[, c(cancer_self_report, cancer_self_report_age, 
                                                cancer_icd10, cancer_icd10_age)]
field_id_cancer <- field_id[c(1:13), ]

# Initialize matrices for case status and age at diagnosis
trait_cancer_out <- matrix(0, nrow = sample_size, ncol = nrow(field_id_cancer))
trait_cancer_age <- matrix(NA, nrow = sample_size, ncol = nrow(field_id_cancer))

for (tt in 1:nrow(field_id_cancer)) {
  # Initialize vectors for this disease
  out_s <- rep(0, sample_size)
  age_s <- rep(NA, sample_size)
  
  # Process self-reported data
  self_dat <- NULL
  for (ss in 1:sample_size) {
    cancer_s <- trait_binary_dat_cancer[ss, c(cancer_self_report, cancer_self_report_age)]
    cnd_s <- any(cancer_s[, cancer_self_report] %in% field_id_cancer$`Self-reported Trait`[tt])
    if (cnd_s) {
      cancer_idx_s <- which(cancer_s[, cancer_self_report] %in% field_id_cancer$`Self-reported Trait`[tt])
      cancer_age_s <- cancer_s[, (cancer_idx_s + length(cancer_self_report))]
      if (length(cancer_idx_s) >= 2) {
        cancer_age_s <- cancer_age_s[1, 1]
      }
      self_dat <- rbind(self_dat, c(ss, cancer_age_s))
    }
  }
  
  # Process ICD10 data
  icd10_dat <- NULL
  for (ss in 1:sample_size) {
    cancer_s <- trait_binary_dat_cancer[ss, c(cancer_icd10, cancer_icd10_age)]
    cnd_s <- any(grepl(field_id_cancer$`ICD10`[tt], cancer_s[, cancer_icd10]))
    if (cnd_s) {
      cancer_idx_s <- which(grepl(field_id_cancer$`ICD10`[tt], cancer_s[, cancer_icd10]))
      cancer_age_s <- cancer_s[, (cancer_idx_s + length(cancer_icd10))]
      if (length(cancer_idx_s) >= 2) {
        cancer_age_s <- min(cancer_age_s)
      }
      icd10_dat <- rbind(icd10_dat, c(ss, cancer_age_s))
    }
  }
  
  # Combine data from both sources
  if (!is.null(icd10_dat) | !is.null(self_dat)) {
    inter_idx <- intersect(self_dat[, 1], icd10_dat[, 1])
    if (length(inter_idx) > 0) {
      event_dat <- rbind(self_dat[!self_dat[, 1] %in% inter_idx, ], icd10_dat)
    } else {
      event_dat <- rbind(self_dat, icd10_dat)
    }
    
    # Mark cases before recruitment
    del_cnd <- ifelse(event_dat[, 2] < recruitment_age[event_dat[, 1]], 1, 0)
    event_dat <- cbind(event_dat, del_cnd)
    
    # Assign case status and age
    out_s[event_dat[, 1]] <- 1
    out_s[event_dat[event_dat[, 3] == 1, 1]] <- NA
    
    age_s[event_dat[, 1]] <- event_dat[, 2]
    age_s[event_dat[event_dat[, 3] == 1, 1]] <- NA
    
    cat("For", field_id_cancer$`Self-reported Trait`[tt], ", ",  
        sum(out_s, na.rm = TRUE), "cases are selected.\n")
  }
  
  trait_cancer_out[, tt] <- out_s
  trait_cancer_age[, tt] <- age_s
}

# Identify non-cancer and extract ages
trait_binary_dat_noncancer <- trait_binary_dat[, c("p41202", "p41204", 
                                                   noncancer_self_report, noncancer_self_report_age)]
field_id_noncancer <- field_id[-c(1:13), ]

# Initialize matrices for case status and age at diagnosis
trait_noncancer_out <- matrix(0, nrow = sample_size, ncol = nrow(field_id_noncancer))
trait_noncancer_age <- matrix(NA, nrow = sample_size, ncol = nrow(field_id_noncancer))

for (tt in 1:nrow(field_id_noncancer)) {
  cat("Processing:", field_id_noncancer$`Self-reported Trait`[tt], "\n")
  
  # Initialize vectors for this disease
  out_s <- rep(0, sample_size)
  age_s <- rep(NA, sample_size)
  
  # Process self-reported data
  self_dat <- NULL
  self_trait <- field_id_noncancer$`Self-reported Trait`[tt]
  if (!is.na(self_trait)) {
    for (ss in 1:sample_size) {
      noncancer_s <- trait_binary_dat_noncancer[ss, c(noncancer_self_report, noncancer_self_report_age)]
      cnd_s <- any(noncancer_s[, noncancer_self_report] %in% self_trait)
      if (cnd_s) {
        noncancer_idx_s <- which(noncancer_s[, noncancer_self_report] %in% self_trait)
        noncancer_age_s <- noncancer_s[, (noncancer_idx_s + length(noncancer_self_report))]
        if (length(noncancer_idx_s) >= 2) {
          noncancer_age_s <- noncancer_age_s[1, 1]
        }
        if (noncancer_age_s != -1) {
          self_dat <- rbind(self_dat, c(ss, noncancer_age_s))
        }
      }
    }
  }
  
  # Process ICD10 data
  icd10_dat <- NULL
  noncancer_date_df_s <- strsplit(field_id_noncancer$`ICD10`[tt], "\\|")[[1]] %>% 
    substr(1, 3) %>% 
    aaply(., 1, function(ss) noncancer_id_date[grep(ss, noncancer_id_date[, 2]), 1]) %>%
    noncancer_date_df[, ., drop = FALSE]
  
  if (ncol(noncancer_date_df_s) != 0) {
    # Main diagnosis
    idx_main <- strsplit(field_id_noncancer$`ICD10`[tt], "\\|")[[1]] %>% 
      alply(., 1, function(ss) grepl(ss, trait_binary_dat_noncancer[, "p41202"])) %>%
      do.call("cbind", .)
    
    main_dat <- NULL
    for (ss in 1:sample_size) {
      if (any(idx_main[ss, ])) {
        date_s <- noncancer_date_df_s[ss, , drop = FALSE]
        date_s2 <- date_s[1, idx_main[ss, ]]
        if (length(date_s2) >= 2) {
          date_s2 <- date_s2[1, 1]
        }
        if (!grepl("Code has event date", as.character(date_s2))) {
          age_s_val <- difftime(as.Date(date_s2), 
                                as.Date(as.character(trait_binary_dat$p34[ss]), "%Y"), 
                                units = "days") %>% as.numeric
          age_s_val <- round(age_s_val/365, 1)
          main_dat <- rbind(main_dat, c(ss, age_s_val))
        }
      }
    }
    
    # Second diagnosis
    idx_second <- strsplit(field_id_noncancer$`ICD10`[tt], "\\|")[[1]] %>% 
      alply(., 1, function(ss) grepl(ss, trait_binary_dat_noncancer[, "p41204"])) %>%
      do.call("cbind", .)
    
    second_dat <- NULL
    for (ss in 1:sample_size) {
      if (any(idx_second[ss, ])) {
        date_s <- noncancer_date_df_s[ss, , drop = FALSE]
        date_s2 <- date_s[1, idx_second[ss, ]]
        if (length(date_s2) >= 2) {
          date_s2 <- date_s2[1, 1]
        }
        if (!grepl("Code has event date", as.character(date_s2))) {
          age_s_val <- difftime(as.Date(date_s2), 
                                as.Date(as.character(trait_binary_dat$p34[ss]), "%Y"), 
                                units = "days") %>% as.numeric
          age_s_val <- round(age_s_val/365, 1)
          second_dat <- rbind(second_dat, c(ss, age_s_val))
        }
      }
    }
    
    # Combine main and second diagnosis
    inter_diag <- intersect(second_dat[, 1], main_dat[, 1])
    icd10_dat <- rbind(main_dat, second_dat[!second_dat[, 1] %in% inter_diag, ])
  }
  
  # Combine data from both sources
  if (!is.null(icd10_dat) | !is.null(self_dat)) {
    inter_idx <- intersect(self_dat[, 1], icd10_dat[, 1])
    if (length(inter_idx) > 0) {
      event_dat <- rbind(self_dat[!self_dat[, 1] %in% inter_idx, ], icd10_dat)
    } else {
      event_dat <- rbind(self_dat, icd10_dat)
    }
    
    # Mark cases before recruitment
    del_cnd <- ifelse(event_dat[, 2] < recruitment_age[event_dat[, 1]], 1, 0)
    event_dat <- cbind(event_dat, del_cnd)
    
    # Assign case status and age
    out_s[event_dat[, 1]] <- 1
    out_s[event_dat[event_dat[, 3] == 1, 1]] <- NA
    
    age_s[event_dat[, 1]] <- event_dat[, 2]
    age_s[event_dat[event_dat[, 3] == 1, 1]] <- NA
    
    cat("For", field_id_noncancer$`Reported Trait`[tt], ", ",  
        sum(out_s, na.rm = TRUE), "cases are selected.\n")
  }
  
  trait_noncancer_out[, tt] <- out_s
  trait_noncancer_age[, tt] <- age_s
}

# Combine cancer and non-cancer data
trait_disease_all_out <- cbind(trait_cancer_out, trait_noncancer_out)
trait_disease_all_age <- cbind(trait_cancer_age, trait_noncancer_age)

# Apply sex-specific filters
## All samples
trait_disease_all_dat <- trait_disease_all_out
trait_disease_all_age_dat <- trait_disease_all_age
for (j in 1:ncol(trait_disease_all_dat)) {
  if (field_id$`Sex Specificity`[j] != "All") {
    trait_disease_all_dat[, j] <- NA
    trait_disease_all_age_dat[, j] <- NA
  }
}

## Female samples
female_idx <- which(all_id %in% female_id)
trait_disease_female_dat <- trait_disease_all_out[female_idx, ]
trait_disease_female_age_dat <- trait_disease_all_age[female_idx, ]
for (j in 1:ncol(trait_disease_female_dat)) {
  if (field_id$`Sex Specificity`[j] == "Male") {
    trait_disease_female_dat[, j] <- NA
    trait_disease_female_age_dat[, j] <- NA
  }
}

## Male samples
male_idx <- which(all_id %in% male_id)
trait_disease_male_dat <- trait_disease_all_out[male_idx, ]
trait_disease_male_age_dat <- trait_disease_all_age[male_idx, ]
for (j in 1:ncol(trait_disease_male_dat)) {
  if (field_id$`Sex Specificity`[j] == "Female") {
    trait_disease_male_dat[, j] <- NA
    trait_disease_male_age_dat[, j] <- NA
  }
}

# Set column names
colnames(trait_disease_all_dat) <- colnames(trait_disease_all_age_dat) <- 
  colnames(trait_disease_female_dat) <- colnames(trait_disease_female_age_dat) <- 
  colnames(trait_disease_male_dat) <- colnames(trait_disease_male_age_dat) <- 
  field_id$PHID

# Save results
## Case status
saveRDS(trait_disease_all_dat, file = paste0(PROJ_PATH, "/03_data/out_pheno/disease_trait_all_dat.rds"))
saveRDS(trait_disease_female_dat, file = paste0(PROJ_PATH, "/03_data/out_pheno/disease_trait_female_dat.rds"))
saveRDS(trait_disease_male_dat, file = paste0(PROJ_PATH, "/03_data/out_pheno/disease_trait_male_dat.rds"))

## Age at diagnosis
saveRDS(trait_disease_all_age_dat, file = paste0(PROJ_PATH, "/03_data/out_pheno/disease_age_all_dat.rds"))
saveRDS(trait_disease_female_age_dat, file = paste0(PROJ_PATH, "/03_data/out_pheno/disease_age_female_dat.rds"))
saveRDS(trait_disease_male_age_dat, file = paste0(PROJ_PATH, "/03_data/out_pheno/disease_age_male_dat.rds"))

# Also save as CSV for easy inspection
write.csv(trait_disease_all_age_dat, file = paste0(PROJ_PATH, "/03_data/out_pheno/disease_age_all_dat.csv"), 
          quote = FALSE, row.names = FALSE)
write.csv(trait_disease_female_age_dat, file = paste0(PROJ_PATH, "/03_data/out_pheno/disease_age_female_dat.csv"), 
          quote = FALSE, row.names = FALSE)
write.csv(trait_disease_male_age_dat, file = paste0(PROJ_PATH, "/03_data/out_pheno/disease_age_male_dat.csv"), 
          quote = FALSE, row.names = FALSE)

# Create summary statistics including age information
trait_summary <- data.frame(
  PHID = field_id$PHID,
  Disease_Name = field_id$`Reported Trait`,
  Sex_Specificity = field_id$`Sex Specificity`,
  All_N = NA, All_Cases = NA, All_Controls = NA, All_MeanAge = NA, All_SDAge = NA,
  Female_N = NA, Female_Cases = NA, Female_Controls = NA, Female_MeanAge = NA, Female_SDAge = NA,
  Male_N = NA, Male_Cases = NA, Male_Controls = NA, Male_MeanAge = NA, Male_SDAge = NA
)

for (i in 1:nrow(field_id)) {
  # All samples
  idx_all <- !is.na(trait_disease_all_dat[, i])
  if (any(idx_all)) {
    trait_summary$All_N[i] <- sum(idx_all)
    trait_summary$All_Cases[i] <- sum(trait_disease_all_dat[idx_all, i], na.rm = TRUE)
    trait_summary$All_Controls[i] <- sum(trait_disease_all_dat[idx_all, i] == 0, na.rm = TRUE)
    
    # Calculate mean age for cases only
    case_idx <- idx_all & trait_disease_all_dat[, i] == 1
    if (any(case_idx)) {
      trait_summary$All_MeanAge[i] <- mean(trait_disease_all_age_dat[case_idx, i], na.rm = TRUE)
      trait_summary$All_SDAge[i] <- sd(trait_disease_all_age_dat[case_idx, i], na.rm = TRUE)
    }
  }
  
  # Female samples
  idx_female <- !is.na(trait_disease_female_dat[, i])
  if (any(idx_female)) {
    trait_summary$Female_N[i] <- sum(idx_female)
    trait_summary$Female_Cases[i] <- sum(trait_disease_female_dat[idx_female, i], na.rm = TRUE)
    trait_summary$Female_Controls[i] <- sum(trait_disease_female_dat[idx_female, i] == 0, na.rm = TRUE)
    
    case_idx_female <- idx_female & trait_disease_female_dat[, i] == 1
    if (any(case_idx_female)) {
      trait_summary$Female_MeanAge[i] <- mean(trait_disease_female_age_dat[case_idx_female, i], na.rm = TRUE)
      trait_summary$Female_SDAge[i] <- sd(trait_disease_female_age_dat[case_idx_female, i], na.rm = TRUE)
    }
  }
  
  # Male samples
  idx_male <- !is.na(trait_disease_male_dat[, i])
  if (any(idx_male)) {
    trait_summary$Male_N[i] <- sum(idx_male)
    trait_summary$Male_Cases[i] <- sum(trait_disease_male_dat[idx_male, i], na.rm = TRUE)
    trait_summary$Male_Controls[i] <- sum(trait_disease_male_dat[idx_male, i] == 0, na.rm = TRUE)
    
    case_idx_male <- idx_male & trait_disease_male_dat[, i] == 1
    if (any(case_idx_male)) {
      trait_summary$Male_MeanAge[i] <- mean(trait_disease_male_age_dat[case_idx_male, i], na.rm = TRUE)
      trait_summary$Male_SDAge[i] <- sd(trait_disease_male_age_dat[case_idx_male, i], na.rm = TRUE)
    }
  }
}

write.csv(trait_summary, file = paste0(PROJ_PATH, "/03_data/out_pheno/disease_trait_age_summary.csv"), 
          quote = FALSE, row.names = FALSE)

cat("Files saved in:", paste0(PROJ_PATH, "/03_data/out_pheno/\n"))