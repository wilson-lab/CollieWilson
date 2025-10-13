## ───────────────────────────── 0) INITIALIZE ────────────────────────────────

# Optional clean slate
rm(list = ls())

# Paths
main_dir <- "/Users/mattcollie/HMS Dropbox/Matt Collie/Connectomics/AOTU_MaleCNS"
plot_dir <- file.path(main_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
setwd(main_dir)

# Packages
suppressPackageStartupMessages({
  library(neuprintr)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(gridExtra)
})

# Set parameters
min_syn <- 5

## ───────────────────────── 1) PULL INPUTS BY REGION   ────────────────────────

# metadata
aotu019_meta <- neuprint_get_meta("AOTU019")
aotu019R_meta <- aotu019_meta[grepl("_R", aotu019_meta$side), ]
aotu019L_meta <- aotu019_meta[grepl("_L", aotu019_meta$side), ]

aotu025_meta <- neuprint_get_meta("AOTU025")
aotu025R_meta <- aotu025_meta[grepl("_R", aotu025_meta$side), ]
aotu025L_meta <- aotu025_meta[grepl("_L", aotu025_meta$side), ]

# inputs
inputs_aotu019R <- flywire_partner_summary2(aotu019R_meta, partners = "in", threshold = min_syn)
inputs_aotu019L <- flywire_partner_summary2(aotu019L_meta, partners = "in", threshold = min_syn)

inputs_aotu025R <- flywire_partner_summary2(aotu025R_meta, partners = "in", threshold = min_syn)
inputs_aotu025L <- flywire_partner_summary2(aotu025L_meta, partners = "in", threshold = min_syn)

## ───────────── 2) SUM BY ROI (collapse _L/_R), AVERAGE L/R ─────────────

strip_lr <- function(x) sub("_[LR]$", "", x)  # only removes final _L or _R

sum_and_avg_by_roi <- function(inputs_R, inputs_L) {
  normR <- inputs_R %>%
    mutate(top_np = strip_lr(top_np)) %>%
    filter(!is.na(top_np), top_np != "")
  
  normL <- inputs_L %>%
    mutate(top_np = strip_lr(top_np)) %>%
    filter(!is.na(top_np), top_np != "")
  
  sR <- normR %>%
    group_by(top_np) %>%
    summarise(sum_weight_R = sum(weight, na.rm = TRUE), .groups = "drop")
  
  sL <- normL %>%
    group_by(top_np) %>%
    summarise(sum_weight_L = sum(weight, na.rm = TRUE), .groups = "drop")
  
  full_join(sR, sL, by = "top_np") %>%
    mutate(avg_sum = rowMeans(across(c(sum_weight_R, sum_weight_L)), na.rm = TRUE)) %>%
    filter(avg_sum >= 100) %>%   # omit any ROI with avg < 100
    select(top_np, avg_sum) %>%
    arrange(desc(avg_sum))
}

aotu019_roi <- sum_and_avg_by_roi(inputs_aotu019R, inputs_aotu019L)
aotu025_roi <- sum_and_avg_by_roi(inputs_aotu025R, inputs_aotu025L)


## ───────────────────── 3) PLOT SIDE-BY-SIDE PIE CHARTS ──────────────────────
make_pie <- function(df, title_txt) {
  df <- df %>%
    mutate(frac = avg_sum / sum(avg_sum), label = top_np)
  
  ggplot(df, aes(x = "", y = avg_sum, fill = top_np)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar(theta = "y") +
    guides(fill = guide_legend(title = "ROI", ncol = 1)) +
    labs(title = title_txt, y = NULL, x = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

p019 <- make_pie(aotu019_roi, "AOTU019: Avg summed inputs by ROI")
p025 <- make_pie(aotu025_roi, "AOTU025: Avg summed inputs by ROI")
pie_combo <- grid.arrange(p019, p025, ncol = 2)

## ────────────────────────── 4) SAVE TO plot_dir ─────────────────────────────
png_file <- file.path(plot_dir, "ROI_input_pies_AOTU019_AOTU025.png")
svg_file <- file.path(plot_dir, "ROI_input_pies_AOTU019_AOTU025.svg")
ggsave(png_file, pie_combo, width = 10, height = 5, dpi = 300)
ggsave(svg_file, pie_combo, width = 10, height = 5)
