########################
## Exposure importance Analysis_4
########################
library(survival)
library(CoxR2)
library(data.table)
library(optparse)
library(rms)

# 1. Input parameters
opt_list <- list(
  make_option("--disease", type = "character", default = "PHD01001",
              help = "INPUT: disease code", metavar = "character"),
  make_option("--imp", type = "integer", default = NULL,
              help = "INPUT: imputation index", metavar = "integer")
)
opt <- parse_args(OptionParser(option_list = opt_list))

# Paths
PROJ <- "/public/home/gw_hychu/swj/"
DATA <- paste0(PROJ, "proteohubProject/02_data/")
SIG_PHM <- "/public/home/gw_hychu/swj/web_out_Cox_all/M-D/sig_PHM_all/"
SIG_PHE <- "/public/home/gw_hychu/swj/web_out_Cox_all/E-D/sig_PHE_all/"
OUT <- paste0(PROJ, "web_out/Attribution_Cox4/")

## 2. Load core data
dis <- readRDS(paste0(DATA, "out_pheno_Cox/disease_trait_all_dat.rds"))
follow <- readRDS(paste0(DATA, "follow_data/follow_data_all.rds"))
cov <- readRDS(paste0(DATA, "cov/cov_all.rds"))
prs <- readRDS(paste0(DATA, "PRS_all.rds"))

if (!is.null(opt$imp)) {
  phm_l <- readRDS(paste0(DATA, "out_pheno_attra/imputed_PHM_list_5.rds"))
  phm <- as.matrix(phm_l[[opt$imp]])
  rm(phm_l); gc()
  
  phe_l <- readRDS(paste0(DATA, "out_pheno_attra/imputed_PHE_list_5.rds"))
  phe <- phe_l[[opt$imp]]
  rm(phe_l); gc()
} else {
  phm <- readRDS(paste0(DATA, "out_pheno_attra/measurement_trait_all.rds"))
  phe <- readRDS(paste0(DATA, "out_pheno_attra/exposure_factors_all.rds"))
}

rm_fixed <- c("PHE04018", "PHE04016")
sub_cols <- c("PHE04001", "PHE04002", "PHE04003")
has_sub <- all(sub_cols %in% names(phe))
rm_cols <- if (has_sub) c(rm_fixed, "PHE04004") else rm_fixed
phe <- phe[, !names(phe) %in% rm_cols, drop = FALSE]

f_phm <- paste0(SIG_PHM, opt$disease, "_all.txt")
v_phm <- if (file.exists(f_phm)) {
  intersect(trimws(as.character(fread(f_phm)[[1]])), colnames(phm))
} else character(0)

f_phe <- paste0(SIG_PHE, opt$disease, "_all.txt")
v_phe <- if (file.exists(f_phe)) {
  intersect(trimws(as.character(fread(f_phe)[[1]])), names(phe))
} else character(0)

v_phe <- setdiff(v_phe, rm_cols)

#Potentially collinear variables
pairs <- list(
  c("PHM01009", "PHM01010"),
  c("PHM01024", "PHM01025"),
  c("PHM01028", "PHM01029"),
  c("PHM01032", "PHM01033"),
  c("PHM01036", "PHM01037"),
  c("PHM01040", "PHM01041"),
  c("PHM01044", "PHM01045"),
  c("PHM01046", "PHM01047"),
  c("PHM01050", "PHM01051"),
  c("PHM01054", "PHM01055")
)

for (p in pairs) {
  if (all(p %in% v_phm)) {
    v_phm <- setdiff(v_phm, p[1])
  }
}

expos <- c(v_phm, v_phe)
if (!length(expos)) stop("No significant exposures found for this disease.")

status <- as.numeric(dis[, opt$disease])
time <- as.numeric(follow[[opt$disease]])

prs_col <- paste0(opt$disease, "_Q")
has_prs <- prs_col %in% names(prs)

if (has_prs) {
  pc_cols <- paste0("p22009_a", 1:18)
  dat <- data.frame(
    time = time,
    status = status,
    Sex = cov$Sex,
    Age = cov$Age,
    PRS_Q = as.numeric(as.character(prs[[prs_col]])),
    cov[, pc_cols, drop = FALSE]
  )
  covars <- c(expos, "Sex", "Age", "PRS_Q", pc_cols)
} else {
  dat <- data.frame(
    time = time,
    status = status,
    Sex = cov$Sex,
    Age = cov$Age
  )
  covars <- c(expos, "Sex", "Age")
}

dat <- cbind(dat, as.data.frame(phm[, v_phm, drop = FALSE]))
dat <- cbind(dat, phe[, v_phe, drop = FALSE])

dat <- na.omit(dat)
dat <- dat[dat$time > 0, ]

n <- nrow(dat)
ev <- sum(dat$status)
if (n < 30 || ev < 5) stop("Sample size or events too low.")

f <- as.formula(paste("Surv(time, status) ~", paste(covars, collapse = " + ")))
fit <- coxph(f, data = dat)

# R2
r2 <- as.numeric(coxr2(fit)$rsq)
cidx <- as.numeric(summary(fit)$concordance[1])

# rms
dd <- datadist(dat)
options(datadist = "dd")
fit_cph <- cph(f, data = dat, x = TRUE, y = TRUE, surv = TRUE)

an <- as.data.frame(anova(fit_cph))
if ("TOTAL" %in% rownames(an)) {
  an <- an[rownames(an) != "TOTAL", , drop = FALSE]
}

chi <- grep("^Chi", names(an), value = TRUE)[1]
pcol <- grep("^P", names(an), value = TRUE)[1]
chisq <- as.numeric(an[[chi]])
pval <- as.numeric(an[[pcol]])
pct <- chisq / sum(chisq) * 100

res <- data.frame(
  Variable = rownames(an),
  Chisq = chisq,
  Chisq_Contribution_Pct = pct,
  P_value = pval,
  Model_CoxR2 = r2,
  Model_C_Index = cidx,
  N = n,
  Events = ev,
  PRS_Included = has_prs,
  stringsAsFactors = FALSE
)

out_dir <- paste0(OUT, opt$disease, "/")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
sfx <- if (!is.null(opt$imp)) paste0("_imp", opt$imp) else ""

fwrite(res, file = paste0(out_dir, "cox4_", opt$disease, sfx, ".tsv"), sep = "\t")