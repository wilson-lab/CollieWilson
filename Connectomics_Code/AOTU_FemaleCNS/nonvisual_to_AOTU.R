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
syn_min       <- 10   # min edge weight to be returned by flywire_partner_summary2
syn_min_total <- 25   # min total per cell_type to keep in summaries

## ────────────────────── 1) FETCH AOTU META & INPUTS ─────────────────────────

# AOTU019 meta (L/R)
meta_019_all <- flytable_meta("AOTU019")
meta_019_L   <- meta_019_all[grepl("left",  meta_019_all$side,  ignore.case = TRUE), ]
meta_019_R   <- meta_019_all[grepl("right", meta_019_all$side,  ignore.case = TRUE), ]

# AOTU025 meta (L/R)
meta_025_all <- flytable_meta("AOTU025")
meta_025_L   <- meta_025_all[grepl("left",  meta_025_all$side,  ignore.case = TRUE), ]
meta_025_R   <- meta_025_all[grepl("right", meta_025_all$side,  ignore.case = TRUE), ]

# Monosynaptic inputs with at least syn_min weight
inputs_019_L <- flywire_partner_summary2(meta_019_L, partners = "in", threshold = syn_min)
inputs_019_R <- flywire_partner_summary2(meta_019_R, partners = "in", threshold = syn_min)
inputs_025_L <- flywire_partner_summary2(meta_025_L, partners = "in", threshold = syn_min)
inputs_025_R <- flywire_partner_summary2(meta_025_R, partners = "in", threshold = syn_min)

## ─────────────────── 2) SUM PER cell_type & AVERAGE L/R ─────────────────────

# Collapse a partner summary to total weight per presynaptic cell_type.
# Filters: drop NA cell_type, drop "LC" types, and enforce total >= syn_min_total.
sum_by_cell_type <- function(df,
                             cell_type_col = "cell_type",
                             weight_col    = "weight",
                             syn_min_total = 0) {
  df %>%
    filter(!is.na(.data[[cell_type_col]])) %>%
    group_by(.data[[cell_type_col]]) %>%
    summarise(total_weight = sum(.data[[weight_col]], na.rm = TRUE), .groups = "drop") %>%
    filter(!grepl("LC", .data[[cell_type_col]], ignore.case = TRUE),
           total_weight >= syn_min_total) %>%
    arrange(desc(total_weight)) %>%
    rename(cell_type = !!cell_type_col)
}

# Combine L/R, keeping the union of cell_types (missing side → 0),
# then compute avg_weight = (L + R) / 2; re-apply global threshold & LC filter.
avg_LR_tables <- function(df_L, df_R, syn_min_total = 0) {
  full_join(
    df_L %>% rename(total_weight_L = total_weight),
    df_R %>% rename(total_weight_R = total_weight),
    by = "cell_type"
  ) %>%
    mutate(
      total_weight_L = replace_na(total_weight_L, 0),
      total_weight_R = replace_na(total_weight_R, 0),
      avg_weight     = (total_weight_L + total_weight_R) / 2
    ) %>%
    filter(avg_weight >= syn_min_total,
           !grepl("LC", cell_type, ignore.case = TRUE)) %>%
    arrange(desc(avg_weight))
}

# Per-side summaries
sum_019_L <- sum_by_cell_type(inputs_019_L, syn_min_total = syn_min_total)
sum_019_R <- sum_by_cell_type(inputs_019_R, syn_min_total = syn_min_total)
sum_025_L <- sum_by_cell_type(inputs_025_L, syn_min_total = syn_min_total)
sum_025_R <- sum_by_cell_type(inputs_025_R, syn_min_total = syn_min_total)

# Averaged L/R for each AOTU
inputs_019_avg_by_type <- avg_LR_tables(sum_019_L, sum_019_R, syn_min_total = syn_min_total)
inputs_025_avg_by_type <- avg_LR_tables(sum_025_L, sum_025_R, syn_min_total = syn_min_total)

## ─────────────────── 3) BUILD RANKED SETS FOR VIGNETTES ─────────────────────

# Shared cell_types: present in both, ranked by average of w019 & w025
shared_tbl <- inner_join(
  inputs_019_avg_by_type %>% select(cell_type, w019 = avg_weight),
  inputs_025_avg_by_type %>% select(cell_type, w025 = avg_weight),
  by = "cell_type"
) %>%
  mutate(combined_avg = (w019 + w025) / 2) %>%
  arrange(desc(combined_avg)) %>%
  slice_head(n = 10)

# Unique to AOTU019 (anti-join vs 025), top 5 by weight
uniq_019_tbl <- anti_join(
  inputs_019_avg_by_type %>% select(cell_type, w = avg_weight),
  inputs_025_avg_by_type %>% select(cell_type),
  by = "cell_type"
) %>%
  arrange(desc(w)) %>%
  slice_head(n = 5)

# Unique to AOTU025 (anti-join vs 019), top 5 by weight
uniq_025_tbl <- anti_join(
  inputs_025_avg_by_type %>% select(cell_type, w = avg_weight),
  inputs_019_avg_by_type %>% select(cell_type),
  by = "cell_type"
) %>%
  arrange(desc(w)) %>%
  slice_head(n = 5)

# Shared width scaling across all three vignettes
scale_width <- function(w, min_w = 0.8, max_w = 6, global_range = NULL) {
  w[is.na(w)] <- 0
  if (is.null(global_range)) global_range <- range(w, na.rm = TRUE)
  if (!is.finite(diff(global_range)) || diff(global_range) == 0) {
    return(rep((min_w + max_w) / 2, length(w)))
  }
  (w - global_range[1]) / diff(global_range) * (max_w - min_w) + min_w
}

global_range <- range(
  c(shared_tbl$w019, shared_tbl$w025, uniq_019_tbl$w, uniq_025_tbl$w),
  na.rm = TRUE
)

## ───────────────────────── 4) PLOT BUILDERS ─────────────────────────────────

# Two-target vignette: shared presynaptic types → AOTU019 and AOTU025
build_vignette_shared <- function(df, global_range, col_019, col_025) {
  if (nrow(df) == 0)
    return(ggplot() + theme_void() +
             annotate("text", x = 0.5, y = 0.5, label = "No shared cell_types", size = 5))
  
  # Nodes: presyn at y=1, targets at y=0
  x_types <- seq(0.05, 0.95, length.out = max(nrow(df), 2))
  nodes <- bind_rows(
    tibble(name = df$cell_type, x = x_types, y = 1.0),
    tibble(name = c("AOTU019", "AOTU025"), x = c(0.25, 0.75), y = 0.0)
  )
  
  # Edges to each target (drop zero weights)
  edges_019 <- tibble(from = df$cell_type, to = "AOTU019", weight = df$w019, col = col_019)
  edges_025 <- tibble(from = df$cell_type, to = "AOTU025", weight = df$w025, col = col_025)
  join_xy <- function(edges, nodes_df) {
    edges %>%
      left_join(nodes_df %>% select(name, x, y), by = c("from" = "name")) %>%
      rename(x1 = x, y1 = y) %>%
      left_join(nodes_df %>% select(name, x, y), by = c("to" = "name")) %>%
      rename(x2 = x, y2 = y)
  }
  edges_019 <- join_xy(edges_019, nodes) %>% filter(weight > 0)
  edges_025 <- join_xy(edges_025, nodes) %>% filter(weight > 0)
  
  edges_019$linewidth <- scale_width(edges_019$weight, global_range = global_range)
  edges_025$linewidth <- scale_width(edges_025$weight, global_range = global_range)
  
  ggplot() +
    geom_segment(
      data = edges_019,
      aes(x = x1, y = y1, xend = x2, yend = y2),
      size = edges_019$linewidth, color = edges_019$col,
      lineend = "round", arrow = arrow(type = "closed", length = unit(5, "pt"))
    ) +
    geom_segment(
      data = edges_025,
      aes(x = x1, y = y1, xend = x2, yend = y2),
      size = edges_025$linewidth, color = edges_025$col,
      lineend = "round", arrow = arrow(type = "closed", length = unit(5, "pt"))
    ) +
    geom_point(data = nodes, aes(x = x, y = y), shape = 21, fill = "white", size = 2.8) +
    geom_text(data = nodes %>% filter(y == 0.0),
              aes(x = x, y = y - 0.06, label = name), vjust = 1, size = 3.1) +
    geom_text(data = nodes %>% filter(y == 1.0),
              aes(x = x, y = y + 0.05, label = name), vjust = 0, size = 2.6, angle = 25) +
    coord_cartesian(xlim = c(-0.05, 1.05), ylim = c(-0.15, 1.15), expand = FALSE) +
    theme_minimal(base_size = 12) +
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(), plot.margin = margin(10, 15, 10, 15)) +
    ggtitle("Shared: top 10 by avg(w019, w025)")
}

# Single-target vignette: unique presynaptic types → one AOTU
build_vignette_unique <- function(df, target_label, global_range, col_target) {
  if (nrow(df) == 0)
    return(ggplot() + theme_void() +
             annotate("text", x = 0.5, y = 0.5,
                      label = paste("No unique cell_types for", target_label), size = 5))
  
  x_types <- seq(0.05, 0.95, length.out = max(nrow(df), 2))
  nodes <- bind_rows(
    tibble(name = df$cell_type, x = x_types, y = 1.0),
    tibble(name = target_label,  x = 0.5,    y = 0.0)
  )
  
  edges <- tibble(from = df$cell_type, to = target_label, weight = df$w) %>%
    left_join(nodes %>% select(name, x, y), by = c("from" = "name")) %>%
    rename(x1 = x, y1 = y) %>%
    left_join(nodes %>% select(name, x, y), by = c("to" = "name")) %>%
    rename(x2 = x, y2 = y) %>%
    filter(weight > 0)
  
  edges$linewidth <- scale_width(edges$weight, global_range = global_range)
  
  ggplot() +
    geom_segment(
      data = edges,
      aes(x = x1, y = y1, xend = x2, yend = y2),
      size = edges$linewidth, color = col_target,
      lineend = "round", arrow = arrow(type = "closed", length = unit(5, "pt"))
    ) +
    geom_point(data = nodes, aes(x = x, y = y), shape = 21, fill = "white", size = 2.8) +
    geom_text(data = nodes %>% filter(y == 0.0),
              aes(x = x, y = y - 0.06, label = name), vjust = 1, size = 3.1) +
    geom_text(data = nodes %>% filter(y == 1.0),
              aes(x = x, y = y + 0.05, label = name), vjust = 0, size = 2.6, angle = 25) +
    coord_cartesian(xlim = c(-0.05, 1.05), ylim = c(-0.15, 1.15), expand = FALSE) +
    theme_minimal(base_size = 12) +
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(), plot.margin = margin(10, 15, 10, 15)) +
    ggtitle(paste0("Unique to ", target_label, ": top 5"))
}

## ───────────────────────── 5) DRAW & SAVE PANELS ────────────────────────────

p_shared <- build_vignette_shared(shared_tbl, global_range, col_019, col_025)
p_u019   <- build_vignette_unique(uniq_019_tbl, "AOTU019", global_range, col_019)
p_u025   <- build_vignette_unique(uniq_025_tbl, "AOTU025", global_range, col_025)

# Width-scaling note (top-right)
max_w <- round(max(global_range, na.rm = TRUE), 1)
annot <- grid::textGrob(
  paste0("Line width scaled to max weight = ", max_w),
  x = 0.98, y = 0.97, just = "right",
  gp = grid::gpar(fontsize = 10, fontface = "italic")
)

final_3 <- grid.arrange(grid.arrange(p_shared, p_u019, p_u025, ncol = 3), top = annot)

# Save SVG & PNG
plotname <- file.path(plot_dir, "sharedinputs_to_aotu")

svglite::svglite(paste0(plotname, ".svg"), width = 22, height = 8)
grid.arrange(grid.arrange(p_shared, p_u019, p_u025, ncol = 3), top = annot)
dev.off()

png(paste0(plotname, ".png"), units = "in", width = 22, height = 8, res = 300)
grid.arrange(grid.arrange(p_shared, p_u019, p_u025, ncol = 3), top = annot)
dev.off()
