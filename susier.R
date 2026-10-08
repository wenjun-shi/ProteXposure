# Load packages
library(data.table)
library(dplyr)
library(susieR)
library(optparse)

# Input parameters
args_list <- list(
  make_option("--pro_g", type = "numeric", default = NULL,
              help = "INPUT: protein group", metavar = "character")
)
opt_parser <- OptionParser(option_list = args_list)
opt <- parse_args(opt_parser)
# Set parameters
PATH <- "/public/home/gw_hychu/swj/geno/web_out/susieR/"
REF_PATH <- "/public/home/gw_hychu/swj/geno/reference/2k_ukb/"
BLOCK_PATH <- "/public/home/gw_hychu/swj/geno/LD_block/genome_block/EUR"

# Load protein list and select group
PRO <- fread("/public/home/gw_hychu/swj/geno/web_out/protein.txt",
             header = FALSE)[[1]] %>% as.character()
start_g <- seq(1, 2919, 3)
end_g <- seq(3, 2919, 3)
pro_sub <- PRO[start_g[opt$pro_g]:end_g[opt$pro_g]]
# Get block number per chromosome
block_num <- sapply(1:22, function(chrom) {
  fread(paste0(BLOCK_PATH, "/chr", chrom, ".bed"), header = FALSE) %>% nrow
})
# Function to process one LD block
process_block <- function(chrom, b, summstats_chr) {
  corr_chr_b_l <- readRDS(paste0(REF_PATH, "LD_matrix/chr", chrom, "/blk_", b, ".rds"))
  common_rs <- intersect(corr_chr_b_l$rs, summstats_chr$rs)
  if (length(common_rs) == 0) return(NULL)
  # Align summary stats to LD matrix order
  idx <- match(common_rs, summstats_chr$rs)
  summstats_chr_b <- summstats_chr[idx, ]
  if (length(common_rs) > 1) {
    ind_chr_b <- corr_chr_b_l$rs %in% common_rs
    R_chr_b <- as.matrix(corr_chr_b_l$corr)[ind_chr_b, ind_chr_b]
    pip_chr_b <- tryCatch({
      susie_rss(bhat = summstats_chr_b$beta,
                shat = summstats_chr_b$se,
                R = R_chr_b,
                n = 33325,
                L = 10)$pip
    }, error = function(e) NULL)
    summstats_chr_b$pip_susie <- if (is.null(pip_chr_b)) NA else pip_chr_b
  } else {
    summstats_chr_b$pip_susie <- NA
  }
  summstats_chr_b
}

# Main loop
for (pro in pro_sub) {
  cat("Start: ", pro, "\n")
  for (sex in c("all", "female", "male")) {
    summstats <- fread(paste0(PATH, pro, "/output/summ_", sex, ".assoc.txt.gz"))
    res_list <- lapply(1:22, function(chrom) {
      cat("chrom: ", chrom, "\n")
      summstats_chr <- summstats %>% filter(chr == chrom)
      block_res <- lapply(seq_len(block_num[chrom]), function(b) {
        process_block(chrom, b, summstats_chr)
      })
      rbindlist(block_res, fill = TRUE)
    })
    res_pip <- rbindlist(res_list, fill = TRUE)
    # Write output
    out_file <- paste0(PATH, pro, "/output/summ_", sex, "2.assoc.txt")
    fwrite(res_pip, file = out_file, sep = "\t")
    system(paste0("gzip -f ", out_file))
  }
}