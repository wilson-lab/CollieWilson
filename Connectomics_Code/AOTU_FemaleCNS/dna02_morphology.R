## ───────────────────────────── 0) INITIALIZE ────────────────────────────────

# Optional clean slate
rm(list = ls())

# Paths
main_dir <- "/Users/mattcollie/HMS Dropbox/Matt Collie/Connectomics/AOTU_FemaleCNS"
plot_dir <- file.path(main_dir, "plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
setwd(main_dir)

# Packages
suppressPackageStartupMessages({
  library(fafbseg)
  library(elmr)
  library(nat)
  library(nat.flybrains)
  library(nat.ggplot)
  library(dplyr)
  library(ggplot2)
})

# Set colors
col_right <- "#ecb150"
col_left <- "#5053ce"

# Get region surfaces
brain_mesh <- FAFB14.surf
lo_mesh_l <- subset(FAFB14NP.surf, "LO_L")
lo_mesh_r <- subset(FAFB14NP.surf, "LO_R")
aotu_mesh_l <- subset(FAFB14NP.surf, "AOTU_L")
aotu_mesh_r <- subset(FAFB14NP.surf, "AOTU_R")
lal_mesh_l <- subset(FAFB14NP.surf, "LAL_L")
lal_mesh_r <- subset(FAFB14NP.surf, "LAL_R")

# View matrix (front view)
full_view <- structure(c(0.999475896, -0.002692447, -0.03225674, 0,
                       -0.004669891, -0.998103201, -0.06138369, 0,
                       -0.032030359, 0.061502226, -0.99759269, 0,
                       0,         0,         0,         1),
                     .Dim = c(4L, 4L))

## ───────────────────────── 1) PLOT   ────────────────────────

# metadata
dna02_meta <- flytable_meta("DNa02")
dna02_right <- dna02_meta[grepl("right", dna02_meta$side), ]
dna02_left <- dna02_meta[grepl("left", dna02_meta$side), ]

# Fetch meshes
right_mesh <- read_cloudvolume_meshes(dna02_right$root_id)
left_mesh <- read_cloudvolume_meshes(dna02_left$root_id)

# ---- Build a single ggplot with layered geom_neuron calls
roi_col = "grey80"
roi_alpha = 0.3

p <- gganat +
  geom_neuron(
    brain_mesh,
    rotation_matrix = full_view,
    cols  = "grey90",
    alpha = 0.3
  )+
  geom_neuron(
    lo_mesh_l,
    rotation_matrix = full_view,
    cols  = roi_col,
    alpha = roi_alpha
  )+
  geom_neuron(
    lo_mesh_r,
    rotation_matrix = full_view,
    cols  = roi_col,
    alpha = roi_alpha
  )+
  geom_neuron(
    aotu_mesh_l,
    rotation_matrix = full_view,
    cols  = roi_col,
    alpha = roi_alpha
  )+
  geom_neuron(
    aotu_mesh_r,
    rotation_matrix = full_view,
    cols  = roi_col,
    alpha = roi_alpha
  )+
  geom_neuron(
    lal_mesh_l,
    rotation_matrix = full_view,
    cols  = roi_col,
    alpha = roi_alpha
  )+
  geom_neuron(
    lal_mesh_r,
    rotation_matrix = full_view,
    cols  = roi_col,
    alpha = roi_alpha
  )+
  geom_neuron(
    right_mesh,
    rotation_matrix = full_view,
    cols  = col_right,
    alpha = 1
  )+
  geom_neuron(
    left_mesh,
    rotation_matrix = full_view,
    cols  = col_left,
    alpha = 1
  )
p

# ---- Save
ggsave("dna02_morphology.png", plot = p, width = 8, height = 8, dpi = 300)
