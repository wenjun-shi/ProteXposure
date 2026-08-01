
library(tidyverse)
library(ggVolcano)

base_dir <- "/public/home/gw_hychu/swj/Logit_result/web_out_female/P2D"      ## change
subdirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)

for (subdir in subdirs) {
  id <- basename(subdir)
  pattern <- paste0("^effect_size_", id, "_female\\.tsv$")        ## all,male,female
  files <- list.files(subdir, pattern = pattern, full.names = TRUE)
  file <- files[1]  
  data <- read_tsv(file, show_col_types = FALSE)
  data <- data %>%
    mutate(regulation = case_when(
      OR > 1 & P_BH < 0.05 ~ "up",
      OR < 1 & P_BH < 0.05 ~ "down",
      TRUE ~ "ns"
    ))
  label_text <- paste0("Up: ", sum(data$regulation == "up"),
                       "\nDown: ", sum(data$regulation == "down"))
  
  p <- gradual_volcano(data, x = "OR", y = "P_BH", label = "Protein",
                       label_number = 20, output = FALSE) +
    scale_color_gradient2(low = "#0072B5", mid = "grey90", high = "#E41A1C",
                          midpoint = 1, guide = "none") +
    theme(axis.line = element_line(linewidth = 1.2, color = "black"),
          axis.ticks = element_line(linewidth = 1.2),
          axis.text = element_text(face = "bold", size = 11),
          axis.title = element_text(size = 15),
          legend.position = "none") +
    coord_cartesian(xlim = c(min(data$OR) * 0.98, max(data$OR) * 1.02)) +
    labs(x = "Odds Ratio") +
    annotate("text", x = -Inf, y = Inf, label = label_text,
             hjust = -0.05, vjust = 1.3, size = 4.5, fontface = "bold", color = "#ce5b4d")
  outfile <- file.path(subdir, paste0(id, "_volcano_plot.png"))
  ggsave(outfile, p, width = 15, height = 12, units = "cm", dpi = 300)
}

