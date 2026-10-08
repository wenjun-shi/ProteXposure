library(bigreadr)
library(plyr)
library(dplyr)

# Parameters
PROJ_PATH      <- "/public/home/gw_hychu/swj/03_trait/proteohubProject"
TRAIT          <- "/public/home/gw_hychu/swj/03_trait/participant.csv.gz"
NONCANCER_DATE <- "/public/home/gw_hychu/swj/03_trait/noncancer_date.rds"
NONCANCER_ID   <- "/public/home/gw_hychu/swj/03_trait/fieldID_noncancer_date.txt"
OUT_DIR        <- paste0(PROJ_PATH, "/03_data/out_pheno")

all_id    <- fread2(paste0(PROJ_PATH, "/sample_id/all_eid"))[, 1]
male_id   <- fread2(paste0(PROJ_PATH, "/sample_id/male_eid"))[, 1]
female_id <- fread2(paste0(PROJ_PATH, "/sample_id/female_eid"))[, 1]

# Field names 
cancer_self_report       <- c(paste0("p20001_i0_a", 0:5), paste0("p20001_i1_a", 0:5),
                              paste0("p20001_i2_a", 0:5), paste0("p20001_i3_a", 0:5))
cancer_self_report_age   <- c(paste0("p20007_i0_a", 0:5), paste0("p20007_i1_a", 0:5),
                              paste0("p20007_i2_a", 0:5), paste0("p20007_i3_a", 0:5))
cancer_icd10             <- paste0("p40006_i", 0:21)
cancer_icd10_age         <- paste0("p40008_i", 0:21)
noncancer_self_report     <- c(paste0("p20002_i0_a", 0:33), paste0("p20002_i1_a", 0:33),
                               paste0("p20002_i2_a", 0:33), paste0("p20002_i3_a", 0:33))
noncancer_self_report_age <- c(paste0("p20009_i0_a", 0:33), paste0("p20009_i1_a", 0:33),
                               paste0("p20009_i2_a", 0:33), paste0("p20009_i3_a", 0:33))

# Load dat
trait_binary_dat <- fread2(TRAIT, select = c("eid", "p21022", "p34", "p41202", "p41204",
                                             cancer_self_report, cancer_self_report_age,
                                             cancer_icd10, cancer_icd10_age,
                                             noncancer_self_report, noncancer_self_report_age))
trait_binary_dat <- trait_binary_dat[match(all_id, trait_binary_dat$eid), ]
sample_size <- nrow(trait_binary_dat)
recruitment_age  <- trait_binary_dat$p21022

noncancer_date_df <- readRDS(NONCANCER_DATE)
noncancer_date_df <- noncancer_date_df[match(all_id, noncancer_date_df$eid), ]
noncancer_id_date <- fread2(NONCANCER_ID, header = F)

field_id <- fread2(paste0(PROJ_PATH, "/02_data/fliedID_disease_trait.txt"))
field_id$ICD10  <- gsub(",", "\\|", field_id$ICD10)


merge_events <- function(self_dat, icd10_dat, recruitment_age) {
  if (is.null(self_dat) && is.null(icd10_dat)) return(NULL)
  if (!is.null(self_dat) && !is.null(icd10_dat)) {
    inter_idx <- intersect(self_dat[, 1], icd10_dat[, 1])
    event_dat <- if (length(inter_idx) > 0)
      rbind(self_dat[!self_dat[, 1] %in% inter_idx, ], icd10_dat)
    else rbind(self_dat, icd10_dat)
  } else {
    event_dat <- if (!is.null(self_dat)) self_dat else icd10_dat
  }
  del_cnd <- ifelse(event_dat[, 2] < recruitment_age[event_dat[, 1]], 1, 0)
  cbind(event_dat, del_cnd)
}

# Assign events to out / age vectors
assign_events <- function(event_dat, sample_size) {
  out_s <- rep(0, sample_size); age_s <- rep(NA, sample_size)
  if (!is.null(event_dat)) {
    out_s[event_dat[, 1]] <- 1
    out_s[event_dat[event_dat[, 3] == 1, 1]] <- NA
    age_s[event_dat[, 1]] <- event_dat[, 2]
    age_s[event_dat[event_dat[, 3] == 1, 1]] <- NA
  }
  list(out = out_s, age = age_s)
}

# Cancer
trait_binary_dat_cancer <- trait_binary_dat[, c(cancer_self_report, cancer_self_report_age,
                                                cancer_icd10, cancer_icd10_age)]
field_id_cancer <- field_id[1:13, ]

trait_cancer_out <- matrix(0,  sample_size, nrow(field_id_cancer))
trait_cancer_age <- matrix(NA, sample_size, nrow(field_id_cancer))

for (tt in 1:nrow(field_id_cancer)) {
  # --- self-report ---
  self_dat <- NULL
  for (ss in 1:sample_size) {
    cs <- trait_binary_dat_cancer[ss, c(cancer_self_report, cancer_self_report_age)]
    cnd_s <- cs[, cancer_self_report] %in% field_id_cancer$`Self-reported Trait`[tt]
    if (any(cnd_s)) {
      idx <- which(cnd_s)
      a   <- cs[, idx + length(cancer_self_report)]
      if (length(idx) >= 2) a <- a[1, 1]
      self_dat <- rbind(self_dat, c(ss, a))
    }
  }
  # --- ICD10 ---
  icd10_dat <- NULL
  for (ss in 1:sample_size) {
    cs <- trait_binary_dat_cancer[ss, c(cancer_icd10, cancer_icd10_age)]
    cnd_s <- grepl(field_id_cancer$ICD10[tt], cs[, cancer_icd10])
    if (any(cnd_s)) {
      idx <- which(cnd_s)
      a   <- cs[, idx + length(cancer_icd10)]
      if (length(idx) >= 2) a <- min(a)
      icd10_dat <- rbind(icd10_dat, c(ss, a))
    }
  }
  res <- assign_events(merge_events(self_dat, icd10_dat, recruitment_age), sample_size)
  trait_cancer_out[, tt] <- res$out
  trait_cancer_age[, tt] <- res$age
}

# Non-cancer
trait_binary_dat_noncancer <- trait_binary_dat[, c("p41202", "p41204",
                                                   noncancer_self_report, noncancer_self_report_age)]
field_id_noncancer <- field_id[-c(1:13), ]

trait_noncancer_out <- matrix(0,  sample_size, nrow(field_id_noncancer))
trait_noncancer_age <- matrix(NA, sample_size, nrow(field_id_noncancer))

for (tt in 1:nrow(field_id_noncancer)) {
  self_dat   <- NULL
  self_trait <- field_id_noncancer$`Self-reported Trait`[tt]
  if (!is.na(self_trait)) {
    for (ss in 1:sample_size) {
      cs <- trait_binary_dat_noncancer[ss, c(noncancer_self_report, noncancer_self_report_age)]
      cnd_s <- cs[, noncancer_self_report] %in% self_trait
      if (any(cnd_s)) {
        idx <- which(cnd_s)
        a   <- cs[, idx + length(noncancer_self_report)]
        if (length(idx) >= 2) a <- a[1, 1]
        if (a != -1) self_dat <- rbind(self_dat, c(ss, a))
      }
    }
  }
  
  # --- ICD10 ---
  icd10_dat <- NULL
  noncancer_date_df_s <- strsplit(field_id_noncancer$ICD10[tt], "\\|")[[1]] %>%
    substr(1, 3) %>%
    aaply(., 1, function(ss) noncancer_id_date[grep(ss, noncancer_id_date[, 2]), 1]) %>%
    noncancer_date_df[, ., drop = FALSE]
  
  if (ncol(noncancer_date_df_s) != 0) {
    # main diagnosis
    idx_main <- strsplit(field_id_noncancer$ICD10[tt], "\\|")[[1]] %>%
      alply(., 1, function(ss) grepl(ss, trait_binary_dat_noncancer[, "p41202"])) %>%
      do.call("cbind", .)
    main_dat <- NULL
    for (ss in 1:sample_size) {
      if (any(idx_main[ss, ])) {
        date_s2 <- noncancer_date_df_s[ss, idx_main[ss, ], drop = FALSE]
        if (length(date_s2) >= 2) date_s2 <- date_s2[1, 1]
        if (!grepl("Code has event date", as.character(date_s2))) {
          age_s_val <- round(as.numeric(difftime(
            as.Date(date_s2),
            as.Date(as.character(trait_binary_dat$p34[ss]), "%Y"),
            units = "days")) / 365, 1)
          main_dat <- rbind(main_dat, c(ss, age_s_val))
        }
      }
    }
    # second diagnosis
    idx_second <- strsplit(field_id_noncancer$ICD10[tt], "\\|")[[1]] %>%
      alply(., 1, function(ss) grepl(ss, trait_binary_dat_noncancer[, "p41204"])) %>%
      do.call("cbind", .)
    second_dat <- NULL
    for (ss in 1:sample_size) {
      if (any(idx_second[ss, ])) {
        date_s2 <- noncancer_date_df_s[ss, idx_second[ss, ], drop = FALSE]
        if (length(date_s2) >= 2) date_s2 <- date_s2[1, 1]
        if (!grepl("Code has event date", as.character(date_s2))) {
          age_s_val <- round(as.numeric(difftime(
            as.Date(date_s2),
            as.Date(as.character(trait_binary_dat$p34[ss]), "%Y"),
            units = "days")) / 365, 1)
          second_dat <- rbind(second_dat, c(ss, age_s_val))
        }
      }
    }
    inter_diag <- intersect(second_dat[, 1], main_dat[, 1])
    icd10_dat  <- rbind(main_dat, second_dat[!second_dat[, 1] %in% inter_diag, ])
  }
  
  res <- assign_events(merge_events(self_dat, icd10_dat, recruitment_age), sample_size)
  trait_noncancer_out[, tt] <- res$out
  trait_noncancer_age[, tt] <- res$age
}

# Combine and sex-specific filtering
trait_disease_all_out <- cbind(trait_cancer_out, trait_noncancer_out)
trait_disease_all_age <- cbind(trait_cancer_age, trait_noncancer_age)

# All samples
trait_disease_all_dat     <- trait_disease_all_out
trait_disease_all_age_dat <- trait_disease_all_age
for (j in 1:ncol(trait_disease_all_dat)) {
  if (field_id$`Sex Specificity`[j] != "All") {
    trait_disease_all_dat[, j]     <- NA
    trait_disease_all_age_dat[, j] <- NA
  }
}
# Female samples
female_idx <- which(all_id %in% female_id)
trait_disease_female_dat     <- trait_disease_all_out[female_idx, ]
trait_disease_female_age_dat <- trait_disease_all_age[female_idx, ]
for (j in 1:ncol(trait_disease_female_dat)) {
  if (field_id$`Sex Specificity`[j] == "Male") {
    trait_disease_female_dat[, j]     <- NA
    trait_disease_female_age_dat[, j] <- NA
  }
}
# Male samples
male_idx <- which(all_id %in% male_id)
trait_disease_male_dat     <- trait_disease_all_out[male_idx, ]
trait_disease_male_age_dat <- trait_disease_all_age[male_idx, ]
for (j in 1:ncol(trait_disease_male_dat)) {
  if (field_id$`Sex Specificity`[j] == "Female") {
    trait_disease_male_dat[, j]     <- NA
    trait_disease_male_age_dat[, j] <- NA
  }
}

# Unify column names
cn <- field_id$PHID
colnames(trait_disease_all_dat) <- colnames(trait_disease_all_age_dat) <- cn
colnames(trait_disease_female_dat) <- colnames(trait_disease_female_age_dat) <- cn
colnames(trait_disease_male_dat) <- colnames(trait_disease_male_age_dat) <- cn

# Save
saveRDS(trait_disease_all_dat,          file.path(OUT_DIR, "disease_trait_all_dat.rds"))
saveRDS(trait_disease_female_dat,       file.path(OUT_DIR, "disease_trait_female_dat.rds"))
saveRDS(trait_disease_male_dat,         file.path(OUT_DIR, "disease_trait_male_dat.rds"))
saveRDS(trait_disease_all_age_dat,      file.path(OUT_DIR, "disease_age_all_dat.rds"))
saveRDS(trait_disease_female_age_dat,   file.path(OUT_DIR, "disease_age_female_dat.rds"))
saveRDS(trait_disease_male_age_dat,     file.path(OUT_DIR, "disease_age_male_dat.rds"))

# Summary statistics
trait_summary <- data.frame(
  PHID = field_id$PHID,
  Disease_Name = field_id$`Reported Trait`,
  Sex_Specificity = field_id$`Sex Specificity`,
  All_N = NA, All_Cases = NA, All_Controls = NA, All_MeanAge = NA, All_SDAge = NA,
  Female_N = NA, Female_Cases = NA, Female_Controls = NA, Female_MeanAge = NA, Female_SDAge = NA,
  Male_N = NA, Male_Cases = NA, Male_Controls = NA, Male_MeanAge = NA, Male_SDAge = NA
)

for (i in 1:nrow(field_id)) {
  # All
  idx_all <- !is.na(trait_disease_all_dat[, i])
  if (any(idx_all)) {
    trait_summary$All_N[i]        <- sum(idx_all)
    trait_summary$All_Cases[i]    <- sum(trait_disease_all_dat[idx_all, i], na.rm = TRUE)
    trait_summary$All_Controls[i] <- sum(trait_disease_all_dat[idx_all, i] == 0, na.rm = TRUE)
    case_idx <- idx_all & trait_disease_all_dat[, i] == 1
    if (any(case_idx)) {
      trait_summary$All_MeanAge[i] <- mean(trait_disease_all_age_dat[case_idx, i], na.rm = TRUE)
      trait_summary$All_SDAge[i]   <- sd(trait_disease_all_age_dat[case_idx, i], na.rm = TRUE)
    }
  }
  # Female
  idx_female <- !is.na(trait_disease_female_dat[, i])
  if (any(idx_female)) {
    trait_summary$Female_N[i]        <- sum(idx_female)
    trait_summary$Female_Cases[i]    <- sum(trait_disease_female_dat[idx_female, i], na.rm = TRUE)
    trait_summary$Female_Controls[i] <- sum(trait_disease_female_dat[idx_female, i] == 0, na.rm = TRUE)
    case_idx_female <- idx_female & trait_disease_female_dat[, i] == 1
    if (any(case_idx_female)) {
      trait_summary$Female_MeanAge[i] <- mean(trait_disease_female_age_dat[case_idx_female, i], na.rm = TRUE)
      trait_summary$Female_SDAge[i]   <- sd(trait_disease_female_age_dat[case_idx_female, i], na.rm = TRUE)
    }
  }
  # Male
  idx_male <- !is.na(trait_disease_male_dat[, i])
  if (any(idx_male)) {
    trait_summary$Male_N[i]        <- sum(idx_male)
    trait_summary$Male_Cases[i]    <- sum(trait_disease_male_dat[idx_male, i], na.rm = TRUE)
    trait_summary$Male_Controls[i] <- sum(trait_disease_male_dat[idx_male, i] == 0, na.rm = TRUE)
    case_idx_male <- idx_male & trait_disease_male_dat[, i] == 1
    if (any(case_idx_male)) {
      trait_summary$Male_MeanAge[i] <- mean(trait_disease_male_age_dat[case_idx_male, i], na.rm = TRUE)
      trait_summary$Male_SDAge[i]   <- sd(trait_disease_male_age_dat[case_idx_male, i], na.rm = TRUE)
    }
  }
}

write.csv(trait_summary, file.path(OUT_DIR, "disease_trait_age_summary.csv"),
          quote = FALSE, row.names = FALSE)