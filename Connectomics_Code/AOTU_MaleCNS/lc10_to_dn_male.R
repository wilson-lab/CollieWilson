###############################################################################
# LC10a → cell-type outputs and cell-type → DN averages (male-cns:v0.9)
# - Fetch LC10a metadata (L/R) and their outputs (thresholded)
# - Summarize outputs by partner (sum_weight, nInputs), drop weak partners
# - For each partner, pull its outputs and record summed weights to target DNs
# - Collapse to per-type L/R means, then average across sides (missing side = 0)
# - Keep top cell types with nonzero DNa02 average
# - Plot: (1) LC10a→type total output; (2) stacked type→DN outputs
#
# CREATED: 10/06/2025 - MC
###############################################################################

## ───────────────────────────── 0) INITIALIZE ────────────────────────────────

# Clear workspace
rm(list=ls())

# Working directory
main_dir <- "/Users/mattcollie/HMS Dropbox/Matt Collie/Connectomics/AOTU_MaleCNS"
setwd(main_dir)

# Output directory for plots
plot_dir <- file.path(main_dir, "plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
data_dir <- file.path(main_dir, "data")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

# Install neuprintr if missing (silent if present)
if (!requireNamespace("neuprintr", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install_github("natverse/neuprintr")
}

# Libraries
suppressPackageStartupMessages({
  library(neuprintr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(tibble)
  library(svglite)
})

# Connectome dataset (male CNS)
neuprintr::neuprint_login(dataset = "male-cns:v0.9")

# Analysis parameters
min_syn     <- 5    # minimum synapse weight to include in LC10a outputs
min_nLC10   <- 25   # minimum number of rows per partner to keep (>= thresholding of LC10a partner strength)
min_synDN   <- 5   # minimum synapse weight to include when checking partner→DN outputs
dn_list     <- c("DNa02","DNa03","DNb02","DNg13","DNa11","DNa01")  # DNs to track


## ────────────────────── 1) FETCH LC10a META & OUTPUTS ───────────────────────

# LC10a metadata (split by hemisphere)
meta_lc10_all <- neuprint_get_meta("LC10a")
meta_lc10_L   <- meta_lc10_all[grepl("_L", meta_lc10_all$name), ]
meta_lc10_R   <- meta_lc10_all[grepl("_R", meta_lc10_all$name), ]

# LC10a outputs (thresholded by min_syn)
outputs_lc10_L <- neuprint_connection_table(meta_lc10_L, partners = "out", threshold = min_syn, details = TRUE)
outputs_lc10_R <- neuprint_connection_table(meta_lc10_R, partners = "out", threshold = min_syn, details = TRUE)


## ───────────────────── 2) HELPERS (SUMMARIZE & FILL DNs) ────────────────────

# Summarize LC10a outputs by partner:
# - For each partner: type (first non-NA), sum_weight (sum of weights), nInputs (row count)
# - Drop partners with nInputs <= min_nLC10
summarize_outputs_by_partner <- function(outputs_tbl, min_nLC10) {
  outputs_tbl %>%
    group_by(partner) %>%
    summarise(
      type       = dplyr::first(type[!is.na(type)]),
      sum_weight = sum(weight, na.rm = TRUE),
      nInputs    = n(),
      .groups    = "drop"
    ) %>%
    filter(nInputs > min_nLC10) %>%
    arrange(desc(sum_weight))
}

# Add DN columns (initialized to 0)
add_dn_columns <- function(df, dn_list) {
  for (dn in dn_list) df[[dn]] <- 0
  df
}

# For each partner row, fetch its OUT connections (thresholded) and record
# summed weights to each DN in dn_list (stores 0 if no match).
# Also count unique output 'type's containing "DN" as total_DN_outputs.
fill_dn_weights_for_summary <- function(summary_df, dn_list, min_synDN) {
  if (nrow(summary_df) == 0) return(summary_df)
  if (!"total_DN_outputs" %in% names(summary_df)) summary_df$total_DN_outputs <- 0L
  
  for (i in seq_len(nrow(summary_df))) {
    pid <- summary_df$partner[i]
    if (is.na(pid)) next
    
    out_tbl <- tryCatch(
      neuprint_connection_table(pid, partners = "out", threshold = min_synDN, details = TRUE),
      error = function(e) NULL
    )
    if (is.null(out_tbl) || nrow(out_tbl) == 0) next
    
    type_weights <- out_tbl %>%
      dplyr::group_by(type) %>%
      dplyr::summarise(w = sum(weight, na.rm = TRUE), .groups = "drop")
    
    # Fill requested DN weights
    for (dn in dn_list) {
      w_dn <- type_weights$w[type_weights$type == dn]
      if (length(w_dn) > 0) summary_df[[dn]][i] <- w_dn[1]
    }
    
    # Count unique output types that contain "DN" (case-insensitive)
    ct <- unique(stats::na.omit(type_weights$type))
    summary_df$total_DN_outputs[i] <- as.integer(sum(grepl("DN", ct, ignore.case = TRUE)))
  }
  summary_df
}

## ───────────────── 3) SUMMARIZE LC10a OUTPUTS & FILL DN WEIGHTS ─────────────

# Partner-level summaries
lc10_L_summary <- summarize_outputs_by_partner(outputs_lc10_L, min_nLC10)
lc10_R_summary <- summarize_outputs_by_partner(outputs_lc10_R, min_nLC10)

# Initialize DN columns
lc10_L_summary <- add_dn_columns(lc10_L_summary, dn_list)
lc10_R_summary <- add_dn_columns(lc10_R_summary, dn_list)

# Fill DN weights (per partner) + total_DN_outputs
lc10_L_summary <- fill_dn_weights_for_summary(lc10_L_summary, dn_list, 15)
lc10_R_summary <- fill_dn_weights_for_summary(lc10_R_summary, dn_list, 15)

## ──────────────── 4) PER-TYPE L/R MEANS & ACROSS-SIDE AVERAGES ──────────────

# --- Left side --------------------------------------------------------------
dn_type_L <- lc10_L_summary %>%
  dplyr::group_by(type) %>%
  dplyr::summarise(
    partner_L  = dplyr::first(partner),
    sum_LC10a_L        = mean(sum_weight, na.rm = TRUE),
    total_DN_outputs_L = mean(total_DN_outputs, na.rm = TRUE),
    dplyr::across(dplyr::all_of(dn_list), ~ mean(.x, na.rm = TRUE), .names = "{.col}_L"),
    .groups = "drop"
  )

# --- Right side -------------------------------------------------------------
dn_type_R <- lc10_R_summary %>%
  dplyr::group_by(type) %>%
  dplyr::summarise(
    partner_R  = dplyr::first(partner),
    sum_LC10a_R        = mean(sum_weight, na.rm = TRUE),
    total_DN_outputs_R = mean(total_DN_outputs, na.rm = TRUE),
    dplyr::across(dplyr::all_of(dn_list), ~ mean(.x, na.rm = TRUE), .names = "{.col}_R"),
    .groups = "drop"
  )

# --- Join L and R, average across sides; keep post id -----------------------
dn_type_avg <- dplyr::full_join(dn_type_L, dn_type_R, by = "type") %>%
  dplyr::mutate(
    partner = dplyr::coalesce(partner_L, partner_R),
    sum_LC10a        = (dplyr::coalesce(sum_LC10a_L, 0)        + dplyr::coalesce(sum_LC10a_R, 0))        / 2,
    total_DN_outputs = (dplyr::coalesce(total_DN_outputs_L, 0) + dplyr::coalesce(total_DN_outputs_R, 0)) / 2
  )

# --- Average DN columns ------------------------------------------------------
for (dn in dn_list) {
  L <- paste0(dn, "_L"); R <- paste0(dn, "_R")
  if (!L %in% names(dn_type_avg)) dn_type_avg[[L]] <- NA_real_
  if (!R %in% names(dn_type_avg)) dn_type_avg[[R]] <- NA_real_
  dn_type_avg[[dn]] <- (dplyr::coalesce(dn_type_avg[[L]], 0) +
                          dplyr::coalesce(dn_type_avg[[R]], 0)) / 2
}

# --- Keep only averaged columns + post id -----------------------------------
dn_type_avg <- dn_type_avg %>%
  dplyr::select(type, partner, sum_LC10a, total_DN_outputs, dplyr::all_of(dn_list))

# --- Filter/order as before -------------------------------------------------
dn_type_avg_DNa02 <- dn_type_avg %>%
  dplyr::filter(DNa02 > 0) %>%
  dplyr::arrange(dplyr::desc(sum_LC10a))

# --- Save output ------------------------------------------------------------
output_file <- file.path(data_dir, "dn_type_avg_DNa02.csv")
write.csv(dn_type_avg_DNa02, output_file, row.names = FALSE)
message("Saved: ", output_file)

## ────────────── 5) PLOT: LC10a TOTAL + STACKED DN OUTPUTS + DN COUNT (TOP 10)

# Top 10 types (keep ordering for display)
top10 <- dn_type_avg_DNa02 %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::mutate(type = factor(type, levels = rev(type)))  # top-down

dn_long <- top10 %>%
  dplyr::select(type, dplyr::all_of(dn_list)) %>%
  tidyr::pivot_longer(cols = dplyr::all_of(dn_list), names_to = "DN", values_to = "weight") %>%
  dplyr::mutate(DN = factor(DN, levels = rev(dn_list)))  # control stacking order

# Round each axis max up to the next multiple of 500
ymax_lc10a <- ceiling(max(top10$sum_LC10a, na.rm = TRUE) / 500) * 500
ymax_dn <- dn_long %>%
  dplyr::group_by(type) %>%
  dplyr::summarise(stack_total = sum(weight, na.rm = TRUE), .groups = "drop") %>%
  dplyr::pull(stack_total) %>%
  max(na.rm = TRUE) %>%
  { ceiling(. / 500) * 500 }

# Round DN-count axis to nearest 5
ymax_dncount <- ceiling(max(top10$total_DN_outputs, na.rm = TRUE) / 5) * 5

# Panel 1: LC10a → cell type
lc10a_plot <- ggplot2::ggplot(top10, ggplot2::aes(x = type, y = sum_LC10a)) +
  ggplot2::geom_col(fill = "grey50") +
  ggplot2::coord_flip() +
  ggplot2::scale_y_continuous(limits = c(0, ymax_lc10a),
                              breaks = seq(0, ymax_lc10a, by = 500)) +
  ggplot2::labs(x = NULL, y = "Avg LC10a output (sum_LC10a)", title = "LC10a → cell type (avg L/R)") +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(hjust = 0.5))

# Panel 2: cell type → DN
dna_stack_plot <- ggplot2::ggplot(dn_long, ggplot2::aes(x = type, y = weight, fill = DN)) +
  ggplot2::geom_col() +
  ggplot2::coord_flip() +
  ggplot2::scale_y_continuous(limits = c(0, ymax_dn),
                              breaks = seq(0, ymax_dn, by = 100)) +
  ggplot2::labs(x = NULL, y = "Avg DN output weight", title = "Cell type → DN (avg L/R)") +
  ggplot2::scale_fill_brewer(palette = "Paired") +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                 legend.position = "right",
                 plot.title = ggplot2::element_text(hjust = 0.5))

# Panel 3: count of DN outputs per type
dncount_plot <- ggplot2::ggplot(top10, ggplot2::aes(x = type, y = total_DN_outputs)) +
  ggplot2::geom_col(fill = "grey50") +
  ggplot2::coord_flip() +
  ggplot2::scale_y_continuous(limits = c(0, ymax_dncount),
                              breaks = seq(0, ymax_dncount, by = 5)) +
  ggplot2::labs(x = NULL, y = "Avg # DN outputs", title = "Count of DN output types (avg L/R)") +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(hjust = 0.5))

# Combine vertically
combined_plot <- gridExtra::grid.arrange(lc10a_plot, dna_stack_plot, dncount_plot, ncol = 1)

## ───────────────────────── 6) SAVE SVG + PNG ────────────────────────────────

plotname <- file.path(plot_dir, "lc10a_to_dn_malesummary")

# SVG
svglite::svglite(paste0(plotname, ".svg"), width = 8, height = 16)
gridExtra::grid.arrange(lc10a_plot, dna_stack_plot, dncount_plot, ncol = 1)
dev.off()

# PNG (high-res)
png(filename = paste0(plotname, ".png"), units = "in", width = 8, height = 16, res = 1200)
gridExtra::grid.arrange(lc10a_plot, dna_stack_plot, dncount_plot, ncol = 1)
dev.off()
