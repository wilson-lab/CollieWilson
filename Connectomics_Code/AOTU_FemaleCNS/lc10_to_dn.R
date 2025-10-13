###############################################################################
# LC10a → cell-type outputs and cell-type → DN averages (female v783)
# - Fetch LC10a metadata (L/R) and their outputs (thresholded)
# - Summarize outputs by partner (sum_weight, nInputs), drop weak partners
# - For each partner, pull its outputs and record summed weights to target DNs
# - Collapse to per-type L/R means, then average across sides (missing side = 0)
# - Keep top cell types with nonzero DNa02 average
# - Plot: (1) LC10a→type total output; (2) stacked type→DN outputs
#
# CREATED: 10/07/2025 - MC
###############################################################################

## ───────────────────────────── 0) INITIALIZE ────────────────────────────────

# Clear workspace
rm(list=ls())

# Working directory
main_dir <- "/Users/mattcollie/HMS Dropbox/Matt Collie/Connectomics/AOTU_FemaleCNS"
setwd(main_dir)

# Output directory for plots
plot_dir <- file.path(main_dir, "plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
data_dir <- file.path(main_dir, "data")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

# Install fafbseg if missing (silent if present)
if (!requireNamespace("fafbseg", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install_github("natverse/fafbseg")
}

# Libraries
suppressPackageStartupMessages({
  library(fafbseg)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(tibble)
  library(svglite)
})

# Connectome dataset (female CNS)
download_flywire_release_data(which = 'all', version = 783L)

# Analysis parameters
min_syn     <- 5    # minimum synapse weight to include in LC10a outputs
min_nLC10   <- 25   # minimum number of rows per partner to keep (>= thresholding of LC10a partner strength)
min_synDN   <- 10   # minimum synapse weight to include when checking partner→DN outputs
dn_list     <- c("DNa02","DNa03","DNb02","DNg13","DNa11","DNa01")  # DNs to track


## ────────────────────── 1) FETCH LC10a META & OUTPUTS ───────────────────────

# LC10a metadata (split by hemisphere)
meta_lc10_all <- flytable_meta("LC10a")
meta_lc10_L   <- meta_lc10_all[grepl("left", meta_lc10_all$side), ]
meta_lc10_R   <- meta_lc10_all[grepl("right", meta_lc10_all$side), ]

# LC10a outputs (thresholded by min_syn)
outputs_lc10_L <- flywire_partner_summary2(meta_lc10_L, partners = "out", threshold = min_syn)
outputs_lc10_R <- flywire_partner_summary2(meta_lc10_R, partners = "out", threshold = min_syn)


## ───────────────────── 2) HELPERS (SUMMARIZE & FILL DNs) ────────────────────

# Summarize LC10a outputs by post_pt_root_id:
# - For each post_pt_root_id: cell_type (first non-NA), sum_weight (sum of weights), nInputs (row count)
# - Drop post_pt_root_id with nInputs <= min_nLC10
summarize_outputs_by_partner <- function(outputs_tbl, min_nLC10) {
  outputs_tbl %>%
    group_by(post_pt_root_id) %>%
    summarise(
      cell_type       = dplyr::first(cell_type[!is.na(cell_type)]),
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

# For each post_pt_root_id row, fetch its OUT connections (thresholded) and record
# summed weights to each DN in dn_list (stores 0 if no match) and the count of
# unique output cell_types containing "DN" as total_DN_outputs.
fill_dn_weights_for_summary <- function(summary_df, dn_list, this_min) {
  if (nrow(summary_df) == 0) return(summary_df)
  if (!"total_DN_outputs" %in% names(summary_df)) summary_df$total_DN_outputs <- 0L
  
  for (i in seq_len(nrow(summary_df))) {
    pid <- summary_df$post_pt_root_id[i]
    if (is.na(pid)) next
    
    out_tbl <- tryCatch(
      flywire_partner_summary2(pid, partners = "out", threshold = this_min),
      error = function(e) NULL
    )
    if (is.null(out_tbl) || nrow(out_tbl) == 0) next
    
    type_weights <- out_tbl %>%
      dplyr::group_by(cell_type) %>%
      dplyr::summarise(w = sum(weight, na.rm = TRUE), .groups = "drop")
    
    # fill requested DN weights
    for (dn in dn_list) {
      w_dn <- type_weights$w[type_weights$cell_type == dn]
      if (length(w_dn) > 0) summary_df[[dn]][i] <- w_dn[1]
    }
    
    # count unique output types that contain "DN"
    ct <- unique(stats::na.omit(type_weights$cell_type))
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

# Fill DN weights (per partner)
lc10_L_summary <- fill_dn_weights_for_summary(lc10_L_summary, dn_list, 15)
lc10_R_summary <- fill_dn_weights_for_summary(lc10_R_summary, dn_list, 15)

## ──────────────── 4) PER-TYPE L/R MEANS & ACROSS-SIDE AVERAGES ──────────────

# --- Left side -------------------------------------------------------------
dn_type_L <- lc10_L_summary %>%
  group_by(cell_type) %>%
  summarise(
    post_pt_root_id_L     = first(post_pt_root_id),   # keep representative ID
    sum_LC10a_L           = mean(sum_weight, na.rm = TRUE),
    total_DN_outputs_L    = mean(total_DN_outputs, na.rm = TRUE),
    across(all_of(dn_list), ~ mean(.x, na.rm = TRUE), .names = "{.col}_L"),
    .groups = "drop"
  )

# --- Right side ------------------------------------------------------------
dn_type_R <- lc10_R_summary %>%
  group_by(cell_type) %>%
  summarise(
    post_pt_root_id_R     = first(post_pt_root_id),
    sum_LC10a_R           = mean(sum_weight, na.rm = TRUE),
    total_DN_outputs_R    = mean(total_DN_outputs, na.rm = TRUE),
    across(all_of(dn_list), ~ mean(.x, na.rm = TRUE), .names = "{.col}_R"),
    .groups = "drop"
  )

# --- Join L and R, average across sides ------------------------------------
dn_type_avg <- full_join(dn_type_L, dn_type_R, by = "cell_type") %>%
  mutate(
    post_pt_root_id  = coalesce(post_pt_root_id_R, post_pt_root_id_L),
    sum_LC10a        = (coalesce(sum_LC10a_L, 0)        + coalesce(sum_LC10a_R, 0))        / 2,
    total_DN_outputs = (coalesce(total_DN_outputs_L, 0) + coalesce(total_DN_outputs_R, 0)) / 2
  )

# --- Average DN columns ----------------------------------------------------
for (dn in dn_list) {
  L <- paste0(dn, "_L"); R <- paste0(dn, "_R")
  if (!L %in% names(dn_type_avg)) dn_type_avg[[L]] <- NA_real_
  if (!R %in% names(dn_type_avg)) dn_type_avg[[R]] <- NA_real_
  dn_type_avg[[dn]] <- (coalesce(dn_type_avg[[L]], 0) + coalesce(dn_type_avg[[R]], 0)) / 2
}

# --- Keep only relevant columns --------------------------------------------
dn_type_avg <- dn_type_avg %>%
  select(cell_type, post_pt_root_id, sum_LC10a, total_DN_outputs, all_of(dn_list))

# --- Restrict to types with nonzero DNa02 and order ------------------------
dn_type_avg_DNa02 <- dn_type_avg %>%
  filter(DNa02 > 0) %>%
  arrange(desc(sum_LC10a))

# --- Save output ------------------------------------------------------------
output_file <- file.path(data_dir, "dn_type_avg_DNa02.csv")
write.csv(dn_type_avg_DNa02, output_file, row.names = FALSE)
message("Saved: ", output_file)

## ────────────── 5) PLOT: LC10a TOTAL + STACKED DN OUTPUTS (TOP 10) ──────────

# Top 10 types (keep ordering for display)
top10 <- dn_type_avg_DNa02 %>%
  slice_head(n = 10) %>%
  mutate(cell_type = factor(cell_type, levels = rev(cell_type)))  # top-down

dn_long <- top10 %>%
  select(cell_type, all_of(dn_list)) %>%
  pivot_longer(cols = all_of(dn_list), names_to = "DN", values_to = "weight") %>%
  mutate(DN = factor(DN, levels = rev(dn_list)))  # control stacking order

# Round each axis max up to the next multiple of 500
ymax_lc10a <- ceiling(max(top10$sum_LC10a, na.rm = TRUE) / 500) * 500
ymax_dn <- dn_long %>%
  dplyr::group_by(cell_type) %>%
  dplyr::summarise(stack_total = sum(weight, na.rm = TRUE), .groups = "drop") %>%
  dplyr::pull(stack_total) %>%
  max(na.rm = TRUE) %>%
  { ceiling(. / 500) * 500 }

# Panel 1: LC10a → cell type
lc10a_plot <- ggplot(top10, aes(x = cell_type, y = sum_LC10a)) +
  geom_col(fill = "grey50") +
  coord_flip() +
  scale_y_continuous(limits = c(0, ymax_lc10a),
                     breaks = seq(0, ymax_lc10a, by = 500)) +
  labs(x = NULL, y = "Avg LC10a output (sum_LC10a)", title = "LC10a → cell type (avg L/R)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(),
        plot.title = element_text(hjust = 0.5))

# Panel 2: cell type → DN
dna_stack_plot <- ggplot(dn_long, aes(x = cell_type, y = weight, fill = DN)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(limits = c(0, ymax_dn),
                     breaks = seq(0, ymax_dn, by = 100)) +
  labs(x = NULL, y = "Avg DN output weight", title = "Cell type → DN (avg L/R)") +
  scale_fill_brewer(palette = "Paired") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5))

# Panel 3: count of DN outputs
# y-limit rounded up to nearest 5
ymax_dncount <- ceiling(max(top10$total_DN_outputs, na.rm = TRUE) / 5) * 5

dncount_plot <- ggplot(top10, aes(x = cell_type, y = total_DN_outputs)) +
  geom_col(fill = "grey50") +
  coord_flip() +
  scale_y_continuous(limits = c(0, ymax_dncount),
                     breaks = seq(0, ymax_dncount, by = 5)) +
  labs(x = NULL, y = "Avg # DN outputs", title = "Count of DN output types (avg L/R)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(),
        plot.title = element_text(hjust = 0.5))

# Combine vertically (3 panels)
combined_plot <- grid.arrange(lc10a_plot, dna_stack_plot, dncount_plot, ncol = 1)


## ───────────────────────── 6) SAVE SVG + PNG ────────────────────────────────

plotname <- file.path(plot_dir, "lc10a_to_dn_femalesummary")

# SVG
svglite::svglite(paste0(plotname, ".svg"), width = 8, height = 16)
gridExtra::grid.arrange(lc10a_plot, dna_stack_plot, dncount_plot, ncol = 1)
dev.off()

# PNG (high-res)
png(filename = paste0(plotname, ".png"), units = "in", width = 8, height = 16, res = 1200)
gridExtra::grid.arrange(lc10a_plot, dna_stack_plot, dncount_plot, ncol = 1)
dev.off()
