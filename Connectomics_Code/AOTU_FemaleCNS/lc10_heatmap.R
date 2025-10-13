###############################################################################
# LC10a → target overview (right hemisphere)
# - For each target neuron type in `noi_list`:
#   * Pull right-side meta, thresholded inputs, keep LC10a presynaptic partners
#   * Fetch LC10a meshes
#   * Plot LC10a meshes over right LO surface (light grey, alpha 0.2)
#   * Save as “[this_neuron]_lc10a_heatmap.png” in plot_dir
#
# CREATED: 10/07/2025 - MC
###############################################################################

## ───────────────────────────── 0) INITIALIZE ────────────────────────────────

# Optional clean slate
rm(list = ls())

# Paths
main_dir <- "/Users/mattcollie/HMS Dropbox/Matt Collie/Connectomics/AOTU_FemaleCNS"
plot_dir <- file.path(main_dir, "heatmaps")
data_dir <- file.path(main_dir, "data")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
setwd(main_dir)

# Packages
suppressPackageStartupMessages({
  library(fafbseg)
  library(nat)
  library(nat.ggplot)
  library(elmr)
  library(nat.flybrains)
  library(dplyr)
  library(ggplot2)
  library(bit64)
})

# Generate heatmap palettes
map2color <- function(x, pal, limits = range(x)){
  max_weight = 60
  pal[findInterval(x, seq(limits[1], max_weight, length.out = length(pal) + 1), 
                   all.inside=TRUE)]
}
ramp_bw = colorRampPalette(c("grey70","black"))
heatmap_bw = ramp_bw(100)
heatmap_col = heat.colors(100)

# Analysis parameters
min_syn  <- 5

# Get region surfaces
lo_mesh <- subset(FAFB14NP.surf, "LO_L")
aotu_mesh <- subset(FAFB14NP.surf, "AOTU_L")

# View matrix (lateral / lobula view)
lo_view <- structure(c(-0.582475244998932, 0.284925609827042, 0.76125580072403, 0,
                       -0.0822701677680016, -0.952387392520905, 0.293520510196686, 0,
                       0.808657646179199, 0.10833965241909, 0.578193068504333, 0,
                       0,         0,         0,         1),
                     .Dim = c(4L, 4L))

# Fetch neurons of interest
dn_type_avg_DNa02 <- read.csv(
  file.path(data_dir, "dn_type_avg_DNa02.csv"),
  colClasses = c(post_pt_root_id = "integer64")
)
message("Loaded: ", file.path(data_dir, "dn_type_avg_DNa02.csv"))

## ───────────────────────── 1) PLOT ALL LC10a CELLS   ────────────────────────

# Right-side metadata for all LC10a neurons
lc10a_meta <- flytable_meta("LC10a")
lc10a_meta <- lc10a_meta[grepl("right", lc10a_meta$side), ]

# Fetch LC10a meshes
all_lc10a_meshes <- read_cloudvolume_meshes(lc10a_meta$root_id)

# ---- Build a single ggplot with layered geom_neuron calls
p <- gganat +
  # Background LO
  geom_neuron(
    lo_mesh,
    rotation_matrix = lo_view,
    cols  = "grey90",
    alpha = 0.3
  ) +
  # Background AOTU
  geom_neuron(
    aotu_mesh,
    rotation_matrix = lo_view,
    cols  = "grey90",
    alpha = 0.3
  ) +
  # LC10a meshes
  geom_neuron(
    all_lc10a_meshes,
    rotation_matrix = lo_view,
    cols = "black",
    alpha = 1
  )
p

# ---- Save
outfile <- file.path(plot_dir,"all_lc10a.png")
ggsave(outfile, plot = p, width = 8, height = 8, dpi = 300)

## ───────────────────────── 2) LOOP & PLOT PER TARGET ────────────────────────

for (n in 1:nrow(dn_type_avg_DNa02)) {
  this_neuron <- dn_type_avg_DNa02$cell_type[n]
  this_ID <- dn_type_avg_DNa02$post_pt_root_id[n]
  message("Processing ", this_neuron, " ...")
  
  # ---- Right-side metadata for the target neuron type
  this_meta <- flytable_meta(this_ID)

  # ---- Thresholded inputs (monosynaptic partners)
  all_inputs <- flywire_partner_summary2(this_meta, partners = "in", threshold = min_syn)
  
  # ---- Keep LC10a inputs only
  lc10a_inputs <- subset(all_inputs, cell_type == "LC10a") %>%
    arrange(weight)
  
  # ---- Fetch LC10a meshes
  lc10a_meshes <- read_cloudvolume_meshes(lc10a_inputs$pre_pt_root_id)
  
  # ---- Build a single ggplot with layered geom_neuron calls (heatmap color)
  p <- gganat +
    # Background LO
    geom_neuron(
      lo_mesh,
      rotation_matrix = lo_view,
      cols  = "grey90",
      alpha = 0.3
    ) +
    # Background AOTU
    geom_neuron(
      aotu_mesh,
      rotation_matrix = lo_view,
      cols  = "grey90",
      alpha = 0.3
    ) +
    # LC10a meshes
    geom_neuron(
      lc10a_meshes,
      rotation_matrix = lo_view,
      cols = map2color(lc10a_inputs$weight, heatmap_col),
      alpha = 1
    )

  # ---- Save
  outfile <- file.path(plot_dir, paste0(this_neuron, "_lc10a_heatmap.png"))
  ggsave(outfile, plot = p, width = 8, height = 8, dpi = 300)
  
  # ---- Build a single ggplot with layered geom_neuron calls (bw)
  p <- gganat +
    # Background LO
    geom_neuron(
      lo_mesh,
      rotation_matrix = lo_view,
      cols  = "white",
      alpha = 1
    ) +
    # Background AOTU
    geom_neuron(
      aotu_mesh,
      rotation_matrix = lo_view,
      cols  = "white",
      alpha = 1
    ) +
    # LC10a meshes
    geom_neuron(
      lc10a_meshes,
      rotation_matrix = lo_view,
      cols = map2color(lc10a_inputs$weight, heatmap_bw),
      alpha = 1
    )
  
  # ---- Save
  outfile <- file.path(plot_dir, paste0(this_neuron, "_lc10a_bw_heatmap.png"))
  ggsave(outfile, plot = p, width = 8, height = 8, dpi = 300)
}
