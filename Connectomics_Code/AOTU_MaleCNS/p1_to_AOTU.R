###############################################################################
# AOTU P1 Connectivity Summary & Directed Plot
# - Fetch AOTU019/025 metadata and inputs (male-cns:v0.9)
# - Count direct P1 inputs
# - Find disynaptic P1→(AOTU inputs) connections
# - Summarize per bodyid and per type, average across L/R (missing side = 0)
# - Draw top-down directed plots (P1 → types → AOTU019/025) with shared scaling
# - Save outputs and print quick console summaries
#
# CREATED: 10/06/2025 - MC
###############################################################################

## ───────────────────────────── 0) INITIALIZE ────────────────────────────────

# Working directory
main_dir <- "/Users/mattcollie/HMS Dropbox/Matt Collie/Connectomics/AOTU_MaleCNS"
setwd(main_dir)

# Output directory (plots)
plot_dir <- file.path(main_dir, "plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Install neuprintr if needed (silent if present)
if (!requireNamespace("neuprintr", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install_github("natverse/neuprintr")
}

# Libraries
suppressPackageStartupMessages({
  library(neuprintr)
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(tibble)
})

# Neuprint dataset login (male central nervous system connectome)
neuprintr::neuprint_login(dataset = "male-cns:v0.9")

# Visual palette
col_019 <- "#75a1e5"  # AOTU019 edges
col_025 <- "#88518f"  # AOTU025 edges
col_p1  <- "#882f2f"  # P1→type edges

# Analysis parameters
syn_min   <- 10  # minimum synapse weight threshold for pulls
plot_thr  <- 15  # minimum weight threshold to draw an edge in the diagram


## ────────────────────── 1) FETCH AOTU META & INPUTS ─────────────────────────

# AOTU019 metadata (split by hemisphere)
meta_019_all <- neuprint_get_meta("AOTU019")
meta_019_L   <- meta_019_all[grepl("L", meta_019_all$name, ignore.case = TRUE), ]
meta_019_R   <- meta_019_all[grepl("R", meta_019_all$name, ignore.case = TRUE), ]

# AOTU025 metadata (split by hemisphere)
meta_025_all <- neuprint_get_meta("AOTU025")
meta_025_L   <- meta_025_all[grepl("L", meta_025_all$name, ignore.case = TRUE), ]
meta_025_R   <- meta_025_all[grepl("R", meta_025_all$name, ignore.case = TRUE), ]

# Inputs to AOTU neurons (monosynaptic partners with at least syn_min weight)
inputs_019_L <- neuprint_connection_table(meta_019_L, partners = "in", threshold = syn_min, details = TRUE)
inputs_019_R <- neuprint_connection_table(meta_019_R, partners = "in", threshold = syn_min, details = TRUE)
inputs_025_L <- neuprint_connection_table(meta_025_L, partners = "in", threshold = syn_min, details = TRUE)
inputs_025_R <- neuprint_connection_table(meta_025_R, partners = "in", threshold = syn_min, details = TRUE)


## ───────────────── 2) QUICK CHECK: DIRECT P1 INPUT COUNTS ───────────────────

count_p1_rows <- function(df) {
  sum(grepl("P1_", df$type, ignore.case = TRUE), na.rm = TRUE)
}

n_019_L_direct_p1 <- count_p1_rows(inputs_019_L)
n_019_R_direct_p1 <- count_p1_rows(inputs_019_R)
n_025_L_direct_p1 <- count_p1_rows(inputs_025_L)
n_025_R_direct_p1 <- count_p1_rows(inputs_025_R)

cat(
  "Direct P1 inputs:\n",
  "  019L (", n_019_L_direct_p1, ")\n",
  "  019R (", n_019_R_direct_p1, ")\n",
  "  025L (", n_025_L_direct_p1, ")\n",
  "  025R (", n_025_R_direct_p1, ")\n",
  sep = ""
)


## ──────────────── 3) DISYNAPTIC: INPUTS-TO-INPUTS (FILTER P1) ───────────────

# Helper: given an inputs_* table, fetch IN partners to its "partner" set
fetch_inputs_to_inputs <- function(inputs_tbl, min_w = syn_min) {
  partner_ids <- unique(na.omit(inputs_tbl$partner))
  if (length(partner_ids) == 0) return(tibble())
  neuprint_connection_table(partner_ids, partners = "in", threshold = min_w, details = TRUE)
}

# Pull inputs→inputs for all four
inputs_inputs_019_L  <- fetch_inputs_to_inputs(inputs_019_L)
inputs_inputs_019_R  <- fetch_inputs_to_inputs(inputs_019_R)
inputs_inputs_025_L  <- fetch_inputs_to_inputs(inputs_025_L)
inputs_inputs_025_R  <- fetch_inputs_to_inputs(inputs_025_R)

# Keep only rows where 'type' matches P1_* (case-insensitive)
inputs_inputs_019_L_P1 <- inputs_inputs_019_L %>% filter(grepl("P1_", type, ignore.case = TRUE))
inputs_inputs_019_R_P1 <- inputs_inputs_019_R %>% filter(grepl("P1_", type, ignore.case = TRUE))
inputs_inputs_025_L_P1 <- inputs_inputs_025_L %>% filter(grepl("P1_", type, ignore.case = TRUE))
inputs_inputs_025_R_P1 <- inputs_inputs_025_R %>% filter(grepl("P1_", type, ignore.case = TRUE))

# Console: unique P1 partners counted per dataset
n_019_L_partners <- length(unique(inputs_inputs_019_L_P1$partner))
n_019_R_partners <- length(unique(inputs_inputs_019_R_P1$partner))
n_025_L_partners <- length(unique(inputs_inputs_025_L_P1$partner))
n_025_R_partners <- length(unique(inputs_inputs_025_R_P1$partner))
cat(
  "Unique P1 partners (disynaptic, by dataset):\n",
  "  019L: ", n_019_L_partners, "\n",
  "  019R: ", n_019_R_partners, "\n",
  "  025L: ", n_025_L_partners, "\n",
  "  025R: ", n_025_R_partners, "\n",
  sep = ""
)


## ───────────── 4) SUMMARIZE P1 INPUT WEIGHT PER TARGET (BY bodyid) ──────────

# For each bodyid in the inputs→inputs P1-filtered tables, sum total P1 weight
summarize_p1_by_bodyid <- function(df_p1) {
  df_p1 %>%
    group_by(bodyid) %>%
    summarise(p1_weight = sum(weight, na.rm = TRUE), .groups = "drop")
}

p1_sum_019_L <- summarize_p1_by_bodyid(inputs_inputs_019_L_P1)
p1_sum_019_R <- summarize_p1_by_bodyid(inputs_inputs_019_R_P1)
p1_sum_025_L <- summarize_p1_by_bodyid(inputs_inputs_025_L_P1)
p1_sum_025_R <- summarize_p1_by_bodyid(inputs_inputs_025_R_P1)


## ───── 5) ENRICH WITH AOTU INPUT METADATA (type, weight at AOTU synapse) ────

# Adds the AOTU input's cell type and AOTU-synapse weight for each bodyid
add_aotu_meta <- function(p1_tbl, aotu_inputs_tbl) {
  meta_tbl <- aotu_inputs_tbl %>%
    select(partner, type, weight) %>%
    distinct() %>%
    mutate(partner = suppressWarnings(as.numeric(partner)))
  
  p1_tbl %>%
    mutate(bodyid = suppressWarnings(as.numeric(bodyid))) %>%
    left_join(meta_tbl, by = c("bodyid" = "partner")) %>%
    relocate(type, weight, .before = 1)  # put meta columns first
}

p1_sum_019_L <- add_aotu_meta(p1_sum_019_L, inputs_019_L)
p1_sum_019_R <- add_aotu_meta(p1_sum_019_R, inputs_019_R)
p1_sum_025_L <- add_aotu_meta(p1_sum_025_L, inputs_025_L)
p1_sum_025_R <- add_aotu_meta(p1_sum_025_R, inputs_025_R)


## ─────── 6) PER-SIDE SUMMARY BY TYPE + L/R-AVERAGED (MISSING = 0) ───────────

# Collapse to per-type per-side means (AOTU weight and P1 weight), then average across sides
per_side_summary <- function(p1_sum_tbl, side = c("L", "R")) {
  side <- match.arg(side)
  if (side == "L") {
    p1_sum_tbl %>%
      group_by(type) %>%
      summarise(weight_L = mean(weight, na.rm = TRUE),
                p1_weight_L = mean(p1_weight, na.rm = TRUE), .groups = "drop")
  } else {
    p1_sum_tbl %>%
      group_by(type) %>%
      summarise(weight_R = mean(weight, na.rm = TRUE),
                p1_weight_R = mean(p1_weight, na.rm = TRUE), .groups = "drop")
  }
}

# AOTU019
t019_L <- per_side_summary(p1_sum_019_L, "L")
t019_R <- per_side_summary(p1_sum_019_R, "R")
summary_019_avg <- full_join(t019_L, t019_R, by = "type") %>%
  mutate(
    weight_mean    = (coalesce(weight_L, 0)    + coalesce(weight_R, 0)) / 2,
    p1_weight_mean = (coalesce(p1_weight_L, 0) + coalesce(p1_weight_R, 0)) / 2
  ) %>%
  select(type,
         weight_L, weight_R, weight_mean,
         p1_weight_L, p1_weight_R, p1_weight_mean)

# AOTU025
t025_L <- per_side_summary(p1_sum_025_L, "L")
t025_R <- per_side_summary(p1_sum_025_R, "R")
summary_025_avg <- full_join(t025_L, t025_R, by = "type") %>%
  mutate(
    weight_mean    = (coalesce(weight_L, 0)    + coalesce(weight_R, 0)) / 2,
    p1_weight_mean = (coalesce(p1_weight_L, 0) + coalesce(p1_weight_R, 0)) / 2
  ) %>%
  select(type,
         weight_L, weight_R, weight_mean,
         p1_weight_L, p1_weight_R, p1_weight_mean)


## ───────────────── 7) PLOTTING HELPERS (DIRECTED GRAPH) ─────────────────────

# Shared width scaling across both panels
scale_width <- function(w, min_w = 0.5, max_w = 6, global_range = NULL) {
  w[is.na(w)] <- 0
  if (is.null(global_range)) global_range <- range(w, na.rm = TRUE)
  if (!is.finite(diff(global_range)) || diff(global_range) == 0) {
    return(rep((min_w + max_w) / 2, length(w)))
  }
  (w - global_range[1]) / diff(global_range) * (max_w - min_w) + min_w
}

# Build one top-down directed plot (P1 → types → target), filtering on thresholds
build_directed_plot <- function(summary_tbl, target_label, global_range, threshold = plot_thr) {
  
  df <- summary_tbl %>%
    mutate(
      weight_mean    = coalesce(weight_mean, 0),
      p1_weight_mean = coalesce(p1_weight_mean, 0)
    ) %>%
    filter(weight_mean > threshold, p1_weight_mean > threshold)
  
  if (nrow(df) == 0) {
    return(
      ggplot() + theme_void() +
        annotate("text", x = 0.5, y = 0.5,
                 label = paste0("No connections > ", threshold, " for ", target_label),
                 size = 5)
    )
  }
  
  # Spread intermediate nodes horizontally
  types   <- df$type
  n_types <- length(types)
  x_types <- seq(0.05, 0.95, length.out = max(n_types, 2))
  
  # Node coordinates
  nodes <- bind_rows(
    tibble(name = "P1",         x = 0.5,  y = 1.0),
    tibble(name = types,        x = x_types, y = 0.5),
    tibble(name = target_label, x = 0.5,  y = 0.0)
  )
  
  # Edges
  edges_in  <- tibble(from = "P1",  to = types,          weight = df$p1_weight_mean)
  edges_out <- tibble(from = types, to = target_label,   weight = df$weight_mean)
  
  # Join coords
  join_xy <- function(edges, nodes_df) {
    edges %>%
      left_join(nodes_df %>% select(name, x, y), by = c("from" = "name")) %>%
      rename(x1 = x, y1 = y) %>%
      left_join(nodes_df %>% select(name, x, y), by = c("to" = "name")) %>%
      rename(x2 = x, y2 = y)
  }
  edges_in  <- join_xy(edges_in,  nodes)
  edges_out <- join_xy(edges_out, nodes)
  
  # Scale widths with shared global range
  edges_in$linewidth  <- scale_width(edges_in$weight,  global_range = global_range)
  edges_out$linewidth <- scale_width(edges_out$weight, global_range = global_range)
  
  # Target-specific color
  col_target <- ifelse(target_label == "AOTU019", col_019, col_025)
  
  ggplot() +
    # P1 → types (P1 color)
    geom_segment(
      data = edges_in,
      aes(x = x1, y = y1, xend = x2, yend = y2),
      color = col_p1, lineend = "round", size = edges_in$linewidth,
      arrow = arrow(type = "closed", length = unit(6, "pt"))
    ) +
    # types → target (AOTU color)
    geom_segment(
      data = edges_out,
      aes(x = x1, y = y1, xend = x2, yend = y2),
      color = col_target, lineend = "round", size = edges_out$linewidth,
      arrow = arrow(type = "closed", length = unit(6, "pt"))
    ) +
    # Nodes and labels
    geom_point(data = nodes, aes(x = x, y = y), shape = 21, fill = "white", size = 3) +
    geom_text(data = nodes, aes(x = x, y = y - 0.04, label = name),
              vjust = 1, size = 2.8) +
    coord_cartesian(xlim = c(-0.1, 1.1), ylim = c(-0.15, 1.1), expand = FALSE) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title  = element_blank(),
      axis.text   = element_blank(),
      axis.ticks  = element_blank(),
      panel.grid  = element_blank(),
      plot.margin = margin(10, 30, 10, 30)
    ) +
    ggtitle(paste0("P1 → Types → ", target_label, " (>", threshold, ")"))
}


## ──────────────── 8) DRAW SIDE-BY-SIDE PLOTS & SAVE PNG/SVG ────────────────

# Shared scaling across both panels (use both AOTUs and both directions)
global_range <- range(
  c(summary_019_avg$weight_mean, summary_025_avg$weight_mean,
    summary_019_avg$p1_weight_mean, summary_025_avg$p1_weight_mean),
  na.rm = TRUE
)

p019 <- build_directed_plot(summary_019_avg, "AOTU019", global_range)
p025 <- build_directed_plot(summary_025_avg, "AOTU025", global_range)

# Add annotation indicating max used for width scaling (top-right)
max_weight <- round(max(global_range, na.rm = TRUE), 1)
annotation <- textGrob(
  paste0("Line width scaled to max weight = ", max_weight),
  x = 0.98, y = 0.97, just = "right",
  gp = gpar(fontsize = 11, fontface = "italic")
)

final_plot <- grid.arrange(grid.arrange(p019, p025, ncol = 2), top = annotation)

# ---- Save as SVG and PNG ----
plotname <- file.path(plot_dir, "p1_to_aotu")

# SVG
svglite::svglite(paste0(plotname, ".svg"), width = 28, height = 10)
grid.arrange(grid.arrange(p019, p025, ncol = 2), top = annotation)
dev.off()

# PNG (high-res)
png(paste0(plotname, ".png"), units = "in", width = 28, height = 10, res = 300)
grid.arrange(grid.arrange(p019, p025, ncol = 2), top = annotation)
dev.off()
