###############################################################################
# AOTU aIPg Connectivity Summary & Directed Plot
# - Count direct aIPg inputs to AOTU019/025 (female v783)
# - Average L/R per cell type
# - Plot as two separate tiles side-by-side (aIPg→AOTU019, aIPg→AOTU025)
# - Shared line-width scaling; no edge when n aIPg inputs == 0
#
# CREATED: 10/06/2025 - MC
###############################################################################

## ───────────────────────────── 0) INITIALIZE ────────────────────────────────

# Clear workspace
rm(list = ls())

# Paths
main_dir <- "/Users/mattcollie/HMS Dropbox/Matt Collie/Connectomics/AOTU_FemaleCNS"
plot_dir <- file.path(main_dir, "plots")
data_dir <- file.path(main_dir, "data")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
setwd(main_dir)

# Libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(svglite)
  library(tibble)
  library(bit64)
  library(fafbseg)
})

# Connectome dataset (female CNS)
download_flywire_release_data(which = "all", version = 783L)

# Visual palette
col_019  <- "#75a1e5"  # AOTU019
col_025  <- "#88518f"  # AOTU025
col_aipg <- "#882f2f"  # aIPg (node color only)

# Analysis parameters
syn_min <- 10   # min synapse weight threshold for partner pulls

## ────────────────────── 1) FETCH AOTU META & INPUTS ─────────────────────────

# AOTU019 (L/R)
meta_019_all <- flytable_meta("AOTU019")
meta_019_L   <- meta_019_all[grepl("left",  meta_019_all$side, ignore.case = TRUE), ]
meta_019_R   <- meta_019_all[grepl("right", meta_019_all$side, ignore.case = TRUE), ]

# AOTU025 (L/R)
meta_025_all <- flytable_meta("AOTU025")
meta_025_L   <- meta_025_all[grepl("left",  meta_025_all$side, ignore.case = TRUE), ]
meta_025_R   <- meta_025_all[grepl("right", meta_025_all$side, ignore.case = TRUE), ]

# Monosynaptic inputs with at least syn_min weight
inputs_019_L <- flywire_partner_summary2(meta_019_L, partners = "in", threshold = syn_min)
inputs_019_R <- flywire_partner_summary2(meta_019_R, partners = "in", threshold = syn_min)
inputs_025_L <- flywire_partner_summary2(meta_025_L, partners = "in", threshold = syn_min)
inputs_025_R <- flywire_partner_summary2(meta_025_R, partners = "in", threshold = syn_min)

# aIPg search hits (Codex), cast IDs as integer64 so %in% works losslessly
search_results <- read.csv(
  file.path(data_dir, "search_results_aipg.csv"),
  colClasses = c(root_id = "integer64")
)
aipg_results <- search_results %>%
  filter(!is.na(cell_type) & grepl("aIPg", as.character(cell_type), ignore.case = FALSE, fixed = TRUE))

search_ids <- aipg_results$root_id

## ─────────────────────── 2) COUNT & SUM aIPg INPUTS ─────────────────────────

# For a partner summary table, sum weight where pre_pt_root_id ∈ aIPg IDs.
# Return (a) total aIPg input weight and (b) count of matching partners.
sum_aipg_inputs <- function(df, pre_col = "pre_pt_root_id", weight_col = "weight") {
  if (!bit64::is.integer64(df[[pre_col]])) {
    df[[pre_col]] <- bit64::as.integer64(df[[pre_col]])
  }
  matches <- dplyr::filter(df, .data[[pre_col]] %in% search_ids)
  tibble(
    `total aIPg input` = if (nrow(matches) > 0) sum(matches[[weight_col]], na.rm = TRUE) else 0,
    `n aIPg inputs`    = nrow(matches)
  )
}

# Per-side summary
aipg_summary <- dplyr::bind_rows(
  `019_L` = sum_aipg_inputs(inputs_019_L),
  `019_R` = sum_aipg_inputs(inputs_019_R),
  `025_L` = sum_aipg_inputs(inputs_025_L),
  `025_R` = sum_aipg_inputs(inputs_025_R),
  .id = "cell"
)

## ─────────────────────── 3) AVERAGE L/R PER CELL TYPE ───────────────────────

# Produce two rows: AOTU019, AOTU025; average L/R for both metrics.
aipg_avg <- aipg_summary %>%
  mutate(
    cell_type = if_else(grepl("019", cell), "AOTU019", "AOTU025"),
    `total aIPg input` = tidyr::replace_na(`total aIPg input`, 0),
    `n aIPg inputs`    = tidyr::replace_na(`n aIPg inputs`, 0)
  ) %>%
  group_by(cell_type) %>%
  summarise(
    `total aIPg input` = mean(`total aIPg input`, na.rm = TRUE),
    `n aIPg inputs`    = mean(`n aIPg inputs`, na.rm = TRUE),
    .groups = "drop"
  )

print(aipg_avg)

## ───────────────────────── 4) PLOTTING HELPERS ──────────────────────────────

# Shared width scaling (handles all-zero safely)
scale_width <- function(w, min_w = 0.5, max_w = 6, global_range = NULL) {
  w[is.na(w)] <- 0
  if (is.null(global_range)) global_range <- range(w, na.rm = TRUE)
  if (!is.finite(diff(global_range)) || diff(global_range) == 0) {
    return(rep((min_w + max_w) / 2, length(w)))
  }
  (w - global_range[1]) / diff(global_range) * (max_w - min_w) + min_w
}

# Single-tile plot builder: aIPg → (target one of AOTU019/025)
# - Draw no edge if n == 0 (requested behavior)
build_aipg_tile <- function(avg_tbl, target_label,
                            global_range = NULL,
                            col_target   = "#75a1e5",
                            col_aipg_node = "#882f2f") {
  
  # Nodes (top-down)
  nodes <- tibble(
    name = c("aIPg", target_label),
    x    = c(0.5, 0.5),
    y    = c(1.0, 0.0)
  )
  
  # Pick the single target row
  row <- avg_tbl %>% filter(cell_type == target_label)
  if (nrow(row) == 0) {
    return(ggplot() + theme_void() +
             annotate("text", x = 0.5, y = 0.5,
                      label = paste("No data for", target_label), size = 5))
  }
  
  # Build edge; suppress if n == 0
  edge_tbl <- tibble(
    from   = "aIPg",
    to     = target_label,
    weight = row$`total aIPg input`[1],
    n_mean = row$`n aIPg inputs`[1]
  ) %>%
    left_join(nodes, by = c("from" = "name")) %>%
    rename(x1 = x, y1 = y) %>%
    left_join(nodes, by = c("to" = "name")) %>%
    rename(x2 = x, y2 = y) %>%
    mutate(linewidth = scale_width(weight, global_range = global_range))
  
  # Start canvas
  p <- ggplot() +
    # Nodes
    geom_point(data = nodes, aes(x = x, y = y), shape = 21, fill = "white", size = 3) +
    geom_text(data = nodes, aes(x = x, y = y - 0.05, label = name),
              vjust = 1, size = 3) +
    coord_cartesian(xlim = c(0, 1), ylim = c(-0.1, 1.1), expand = FALSE) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(10, 30, 10, 30)
    ) +
    ggtitle(paste0("aIPg → ", target_label, " (avg L/R)"))
  
  # Edge only if n_mean > 0
  if (is.finite(edge_tbl$n_mean) && edge_tbl$n_mean > 0) {
    p <- p +
      geom_segment(
        data = edge_tbl,
        aes(x = x1, y = y1, xend = x2, yend = y2),
        size = edge_tbl$linewidth, color = col_target,
        lineend = "round",
        arrow = arrow(type = "closed", length = unit(6, "pt"))
      ) +
      geom_text(
        data = edge_tbl,
        aes(x = (x1 + x2) / 2, y = (y1 + y2) / 2 + 0.06,
            label = paste0("n=", round(n_mean, 1))),
        size = 3
      )
  } else {
    p <- p +
      annotate("text", x = 0.5, y = 0.5, label = "n=0", size = 3.2, fontface = "italic")
  }
  
  p
}

## ───────────────────────── 5) BUILD TILES & SAVE ────────────────────────────

# Shared scaling across both targets
global_range <- range(aipg_avg$`total aIPg input`, na.rm = TRUE)

# Two separate tiles (side-by-side)
p019 <- build_aipg_tile(aipg_avg, "AOTU019", global_range, col_target = col_019)
p025 <- build_aipg_tile(aipg_avg, "AOTU025", global_range, col_target = col_025)

# Add a small annotation about width scaling
max_w <- round(max(global_range, na.rm = TRUE), 1)
annot <- grid::textGrob(
  paste0("Line width scaled to max weight = ", max_w),
  x = 0.98, y = 0.97, just = "right",
  gp = grid::gpar(fontsize = 11, fontface = "italic")
)

# Arrange tiles; save SVG and PNG
final_plot <- grid.arrange(grid.arrange(p019, p025, ncol = 2), top = annot)

plotname <- file.path(plot_dir, "aipg_to_aotu")

svglite::svglite(paste0(plotname, ".svg"), width = 16, height = 8)
grid.arrange(grid.arrange(p019, p025, ncol = 2), top = annot)
dev.off()

png(paste0(plotname, ".png"), units = "in", width = 16, height = 8, res = 300)
grid.arrange(grid.arrange(p019, p025, ncol = 2), top = annot)
dev.off()
