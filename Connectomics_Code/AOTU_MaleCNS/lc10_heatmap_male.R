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
#main_dir <- "C:/Users/wilson/HMS Dropbox/Matt Collie/Connectomics/AOTU_MaleCNS"
main_dir <- "/Users/mattcollie/HMS Dropbox/Matt Collie/Connectomics/AOTU_MaleCNS"
plot_dir <- file.path(main_dir, "heatmaps")
data_dir <- file.path(main_dir, "data")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
setwd(main_dir)

# Packages
suppressPackageStartupMessages({
  library(neuprintr)
  library(nat.ggplot)
  library(malecns)
  library(dplyr)
  library(ggplot2)
})

# Select male CNS
neuprintr::neuprint_login(dataset = "male-cns:v0.9")

# Generate heatmap palettes
map2color <- function(x, pal, limits = range(x)){
  max_weight = 170
  pal[findInterval(x, seq(limits[1], max_weight, length.out = length(pal) + 1), 
                   all.inside=TRUE)]
}
ramp_bw = colorRampPalette(c("grey70","black"))
heatmap_bw = ramp_bw(100)
heatmap_col = heat.colors(100)

# Analysis parameters
min_syn  <- 5

# Get region surfaces
lo_mesh <- neuprint_ROI_mesh('LO(L)',dataset = "male-cns:v0.9")
aotu_mesh <- neuprint_ROI_mesh('AOTU(L)',dataset = "male-cns:v0.9")

# View matrix (lateral / lobula view)
lo_view <- structure(c(-0.871912717819214, 0.156003504991531, 0.46414577960968, 0,
                       -0.192713916301727, -0.980719804763794, -0.0323898419737816, 0,
                       0.45014414191246, -0.117688581347466, 0.885166466236115, 0,
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
lc10a_meta <- neuprint_get_meta("LC10a")
lc10a_meta <- lc10a_meta[grepl("_L", lc10a_meta$name), ]

# Fetch LC10a meshes
all_lc10a_meshes <- neuprint_read_skeletons(lc10a_meta$bodyid)

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
    root = 0,
    cols = "black",
    alpha = 1
  )

# ---- Save
outfile <- file.path(plot_dir,"all_lc10a.png")
ggsave(outfile, plot = p, width = 8, height = 8, dpi = 300)

## ───────────────────────── 2) LOOP & PLOT PER TARGET ────────────────────────

for (n in 1:nrow(dn_type_avg_DNa02)) {
  this_neuron <- dn_type_avg_DNa02$type[n]
  this_ID <- dn_type_avg_DNa02$partner[n]
  message("Processing ", this_neuron, " ...")

  # ---- Right-side metadata for the target neuron type
  this_meta <- neuprint_get_meta(this_ID)

  # ---- Thresholded inputs (monosynaptic partners)
  all_inputs <- neuprint_connection_table(this_meta, partners = "in", threshold = min_syn, details = TRUE)
  
  # ---- Keep LC10a inputs only
  lc10a_inputs <- subset(all_inputs, type == "LC10a") %>%
    arrange(weight)
  
  # ---- Fetch LC10a meshes
  lc10a_meshes <- neuprint_read_skeletons(lc10a_inputs$partner)
  
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
      root = 0,
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
      root = 0,
      cols = map2color(lc10a_inputs$weight, heatmap_bw),
      alpha = 1
    )
  
  # ---- Save
  outfile <- file.path(plot_dir, paste0(this_neuron, "_lc10a_bw_heatmap.png"))
  ggsave(outfile, plot = p, width = 8, height = 8, dpi = 300)
}
