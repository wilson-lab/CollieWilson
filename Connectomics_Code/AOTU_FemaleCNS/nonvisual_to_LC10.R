###############################################################################
# AOTU inputs by presynaptic cell_type — three vignettes
# - Collapse inputs per presynaptic cell_type (per side), filter, then average L/R
# - Vignette 1: Shared (top 10 by avg of w019 & w025)
# - Vignette 2: Unique to AOTU019 (top 5)
# - Vignette 3: Unique to AOTU025 (top 5)
#
# CREATED: 10/06/2025 - MC
###############################################################################

## ───────────────────────────── 0) INITIALIZE ────────────────────────────────

# Optional clean slate
rm(list = ls())

# Paths
main_dir <- "/Users/mattcollie/HMS Dropbox/Matt Collie/Connectomics/AOTU_FemaleCNS"
plot_dir <- file.path(main_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
setwd(main_dir)

# Libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(svglite)
  library(bit64)
  library(fafbseg)
})

# Connectome dataset (female CNS)
download_flywire_release_data(which = "all", version = 783L)

# Colors
col_019 <- "#75a1e5"  # AOTU019
col_025 <- "#88518f"  # AOTU025

# Analysis parameters
syn_min <- 5   # min edge weight to be returned by flywire_partner_summary2
syn_min_total <- 25

## ────────────────────── 1) FETCH AOTU META & INPUTS ─────────────────────────

# LC10 meta (L/R)
meta_lc10_all <- flytable_meta("LC10a")
meta_lc10_L   <- meta_lc10_all[grepl("left",  meta_lc10_all$side,  ignore.case = TRUE), ]
meta_lc10_R   <- meta_lc10_all[grepl("right", meta_lc10_all$side,  ignore.case = TRUE), ]

# Monosynaptic inputs with at least syn_min weight
inputs_lc10_L <- flywire_partner_summary2(meta_lc10_L, partners = "in", threshold = syn_min)
inputs_lc10_R <- flywire_partner_summary2(meta_lc10_R, partners = "in", threshold = syn_min)

## ─────────────────── 2) SUM PER cell_type & AVERAGE L/R ─────────────────────

# Collapse a partner summary to total weight per presynaptic cell_type.
# Filters: drop NA cell_type, enforce total >= syn_min_total.
sum_by_cell_type <- function(df,
                             cell_type_col = "cell_type",
                             weight_col    = "weight",
                             id_col        = "pre_pt_root_id",
                             syn_min_total = 0) {
  # Ensure ID column is integer64 (safe for n_distinct)
  if (!bit64::is.integer64(df[[id_col]])) {
    df[[id_col]] <- bit64::as.integer64(df[[id_col]])
  }
  
  df %>%
    filter(!is.na(.data[[cell_type_col]])) %>%
    group_by(.data[[cell_type_col]]) %>%
    summarise(
      total_weight  = sum(.data[[weight_col]], na.rm = TRUE),
      n_rows        = dplyr::n(),                                     # (1)
      n_unique_ids  = dplyr::n_distinct(.data[[id_col]]),             # (2)
      .groups = "drop"
    ) %>%
    mutate(
      approx_outputs_per_cell_type = n_rows / pmax(n_unique_ids, 1)   # (3) safe divide
    ) %>%
    filter(total_weight >= syn_min_total) %>%
    arrange(desc(total_weight)) %>%
    rename(cell_type = !!cell_type_col)
}

# Combine L/R, keeping the union of cell_types (missing side → 0),
# then compute avg_weight = (L + R) / 2; re-apply global threshold.
avg_LR_tables <- function(df_L, df_R, syn_min_total = 0) {
  full_join(df_L, df_R, by = "cell_type", suffix = c("_L", "_R")) %>%
    mutate(
      across(
        c(total_weight_L, total_weight_R,
          n_rows_L, n_rows_R,
          n_unique_ids_L, n_unique_ids_R,
          approx_outputs_per_cell_type_L, approx_outputs_per_cell_type_R),
        ~ replace_na(., 0)
      ),
      avg_weight = (total_weight_L + total_weight_R) / 2,
      avg_n_rows = (n_rows_L + n_rows_R) / 2,
      avg_n_ids  = (n_unique_ids_L + n_unique_ids_R) / 2,
      avg_outputs_per_cell_type =
        (approx_outputs_per_cell_type_L + approx_outputs_per_cell_type_R) / 2
    ) %>%
    transmute(
      cell_type,
      avg_weight,
      avg_n_rows,
      avg_n_ids,
      avg_outputs_per_cell_type
    ) %>%
    filter(avg_weight >= syn_min_total) %>%
    arrange(desc(avg_weight))
}


# Per-side summaries
sum_lc10_L <- sum_by_cell_type(inputs_lc10_L, syn_min_total = syn_min_total)
sum_lc10_R <- sum_by_cell_type(inputs_lc10_R, syn_min_total = syn_min_total)

# Averaged L/R
# Compute average per cell type, add weight_per_cell, then filter out LC/Tm
inputs_avg_by_type <- avg_LR_tables(sum_lc10_L, sum_lc10_R, syn_min_total = syn_min_total) %>%
  mutate(
    # Divide each avg_weight by total number of LC10 cells in the meta table
    weight_per_cell = avg_weight / nrow(meta_lc10_all)
  ) %>%
  # Remove rows containing "LC" or "Tm" in the cell_type name
  filter(!grepl("LC", cell_type, ignore.case = TRUE),
         !grepl("Tm", cell_type, ignore.case = TRUE)) %>%
  arrange(desc(weight_per_cell))
