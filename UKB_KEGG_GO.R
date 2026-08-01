#!/bin/bash

########################
## Enrichment for the significant proteins associated with disease
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(plyr)
library(ggplot2)
library(ggprism)
library(gground)

N_PATHWAY <- 5
ALPHA <- 0.05

INPUT_FILE <- "/public/home/gw_hychu/swj/proteohubProject/02_data/sig_proteins_all/sig_proteins_all.tsv"
OUTPUT_DIR <- "/public/home/gw_hychu/swj/web_out_all/enrichment"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)


INPUT_FILE <- "C:/Users/施文俊/Desktop/sig_proteins_all_Logit.tsv"  

INPUT_FILE <- "C:/Users/施文俊/Desktop/sig_proteins_all.tsv"      ##修改
OUTPUT_DIR <- "C:/Users/施文俊/Desktop/web_out_all/enrichment"        ##修改
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

ora_run <- function(pro_symbol){
  gene_symbol <- toupper(pro_symbol)
  gene_id <- bitr(gene_symbol, fromType = "SYMBOL",
                  toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
  
  ora_res_GO <- enrichGO(gene = gene_id[, 1], keyType = 'SYMBOL',
                         OrgDb = org.Hs.eg.db, ont = "ALL",
                         pAdjustMethod = "BH", pvalueCutoff = 1,
                         readable = TRUE)@result
  
  ora_res_KEGG <- enrichKEGG(gene = gene_id[, 2], organism = "hsa",
                             keyType = "ncbi-geneid",
                             pvalueCutoff = 1)@result
  
  ora_res_KEGG$geneID <- aaply(ora_res_KEGG$geneID, 1, function(ss){
    gene_eid <- strsplit(ss, "\\/") %>% unlist 
    gene_sym <- gene_id[gene_id$ENTREZID %in% gene_eid, 1, drop = T] %>% 
      paste0(collapse ="/")
    return(gene_sym)
  })
  
  return(list(GO = ora_res_GO, KEGG = ora_res_KEGG))
}

ora_plt <- function(ora_res){

  go_sig <- ora_res[['GO']] %>% filter(p.adjust < ALPHA)
  kegg_sig <- ora_res[['KEGG']] %>% filter(p.adjust < ALPHA)
  if(nrow(go_sig) == 0 && nrow(kegg_sig) == 0) return(NULL)
  
  use_pathway <- data.frame()
  if(nrow(go_sig) > 0) {
    go_selected <- group_by(go_sig, ONTOLOGY) %>%
      top_n(N_PATHWAY, wt = -p.adjust) %>%
      group_by(p.adjust) %>% top_n(1, wt = Count)
    use_pathway <- rbind(use_pathway, go_selected)
  }
  
  if(nrow(kegg_sig) > 0) {
    kegg_selected <- top_n(kegg_sig, N_PATHWAY, -p.adjust) %>%
      group_by(p.adjust) %>% top_n(1, wt = Count) %>% mutate(ONTOLOGY = 'KEGG')
    use_pathway <- rbind(use_pathway, kegg_selected)
  }
  
  if(nrow(use_pathway) == 0) return(NULL)
  
  use_pathway <- use_pathway %>%
    ungroup() %>%
    mutate(ONTOLOGY = factor(ONTOLOGY, levels = rev(c('BP', 'CC', 'MF', 'KEGG')))) %>%
    dplyr::arrange(ONTOLOGY, p.adjust) %>%
    mutate(Description = factor(Description, levels = Description)) %>%
    tibble::rowid_to_column('index')
  
  gene_out <- strsplit(use_pathway$geneID, "\\/") %>% aaply(., 1, function(x){
    if(length(x) > 6) sample(x, 6) %>% paste0(collapse = "/") %>% paste0("/...(", length(x) -6, ")") 
    else paste0(x, collapse = "/")
  })
  
  width <- 0.5
  xaxis_max <- max(-log10(use_pathway$p.adjust)) + 1
  
  rect.data <- group_by(use_pathway, ONTOLOGY) %>% reframe(n = n()) %>% ungroup() %>%
    mutate(xmin = -3 * width, xmax = -2 * width,
           ymax = cumsum(n), ymin = lag(ymax, default = 0) + 0.6, ymax = ymax + 0.4)
  
  pal <- c('#A0D2E3', '#E79D8A', '#ABD4BA', '#D89FB3')
  
  plt <- ggplot(use_pathway, aes(-log10(p.adjust), y = index, fill = ONTOLOGY)) +
    geom_col(aes(y = Description), width = 0.6, alpha = 0.8) +
    geom_text(aes(x = 0.05, label = Description), hjust = 0, size = 5) +
    geom_text(aes(x = 0.1, label = gene_out, colour = ONTOLOGY), 
              hjust = 0, vjust = 3, size = 3.5, fontface = 'italic', show.legend = FALSE) +
    geom_point(aes(x = -width, size = Count), shape = 21) +
    geom_text(aes(x = -width, label = Count)) +
    scale_size_continuous(name = 'Count', range = c(2, 8)) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = ONTOLOGY),
              data = rect.data, inherit.aes = FALSE) +
    geom_text(aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = ONTOLOGY),
              data = rect.data, angle = 90, hjust = 0.5, vjust = 0.5, inherit.aes = FALSE) +
    annotate("segment", x = 0, y = 0, xend = xaxis_max, yend = 0, linewidth = 1.5) +
    labs(y = NULL, x = expression("-log"[10]*"FDR")) +
    scale_fill_manual(name = 'Category', values = pal) +
    scale_colour_manual(values = pal) +
    scale_x_continuous(breaks = seq(0, xaxis_max, 2), expand = expansion(c(0, 0))) +
    theme_prism() +
    theme(axis.text.y = element_blank(), axis.line = element_blank(),
          axis.ticks.y = element_blank(), legend.title = element_text())
  return(plt)
}


# ORA analysis
data <- read.delim(INPUT_FILE, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

for(i in 1:ncol(data)) {
  disease_name <- colnames(data)[i]
  cat(sprintf("analysis: %s\n", disease_name))
  
  proteins <- data[, i]
  proteins <- proteins[!is.na(proteins) & proteins != ""]
  proteins <- unique(proteins)
  
  if(is.na(proteins[1]) || length(proteins) < 6) {
    # GO
    go_header <- c("ONTOLOGY", "ID", "Description", "GeneRatio", "BgRatio", 
                   "RichFactor", "FoldEnrichment", "zScore",
                   "pvalue", "p.adjust", "qvalue", "geneID", "Count")
    write.table(t(go_header), paste0(OUTPUT_DIR, "/", disease_name, "_enrichment_GO.tsv"),
                sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
    #KEGG
    kegg_header <- c("category", "subcategory", "ID", "Description", "GeneRatio",  
                     "BgRatio", "RichFactor", "FoldEnrichment", "zScore",
                     "pvalue", "p.adjust", "qvalue", "geneID", "Count")
    write.table(t(kegg_header), paste0(OUTPUT_DIR, "/", disease_name, "_enrichment_KEGG.tsv"),
                sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
    next
  }
  
  tryCatch({
    ora_result <- ora_run(proteins)
    # GO
    go_sig <- ora_result$GO %>% filter(p.adjust < ALPHA)
    go_file <- paste0(OUTPUT_DIR, "/", disease_name, "_enrichment_GO.tsv")
    if(nrow(go_sig) > 0) {
      write.table(go_sig, go_file, sep = "\t", quote = FALSE, row.names = FALSE)
    } else {
      go_header <- c("ONTOLOGY", "ID", "Description", "GeneRatio", "BgRatio", 
                     "RichFactor", "FoldEnrichment", "zScore",
                     "pvalue", "p.adjust", "qvalue", "geneID", "Count")
      write.table(t(go_header), go_file, sep = "\t", quote = FALSE,
                  col.names = FALSE, row.names = FALSE)
    }
    # KEGG
    kegg_sig <- ora_result$KEGG %>% filter(p.adjust < ALPHA)
    kegg_file <- paste0(OUTPUT_DIR, "/", disease_name, "_enrichment_KEGG.tsv")
    if(nrow(kegg_sig) > 0) {
      write.table(kegg_sig, kegg_file, sep = "\t", quote = FALSE, row.names = FALSE)
    } else {
      kegg_header <- c("category", "subcategory", "ID", "Description", "GeneRatio",  
                       "BgRatio", "RichFactor", "FoldEnrichment", "zScore",
                       "pvalue", "p.adjust", "qvalue", "geneID", "Count")
      write.table(t(kegg_header), kegg_file, sep = "\t", quote = FALSE,
                  col.names = FALSE, row.names = FALSE)
    }
    
    # png
    if(nrow(go_sig) > 0 || nrow(kegg_sig) > 0) {
      enrich_plot <- ora_plt(ora_result)
      if(!is.null(enrich_plot)) {
        plot_file <- paste0(OUTPUT_DIR, "/", disease_name, "_enrichment_plot.png")
        ggsave(plot_file, enrich_plot, width = 18, height = 12, units = "in", dpi = 320)
      }
    }
    
  }, error = function(e) {
    cat("analysis error\n")
  })
}
warnings()
