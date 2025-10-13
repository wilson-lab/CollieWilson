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
col_AOTU019 <- "#75a1e5"
col_AOTU025 <- "#88518f"

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
aotu019_meta <- flytable_meta("AOTU019")
aotu019_meta <- aotu019_meta[grepl("right", aotu019_meta$side), ]

aotu025_meta <- flytable_meta("AOTU025")
aotu025_meta <- aotu025_meta[grepl("left", aotu025_meta$side), ]

# Fetch meshes
aotu019_mesh <- read_cloudvolume_meshes(aotu019_meta$root_id)
aotu025_mesh <- read_cloudvolume_meshes(aotu025_meta$root_id)

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
    aotu019_mesh,
    rotation_matrix = full_view,
    cols  = col_AOTU019,
    alpha = 1
  )+
  geom_neuron(
    aotu025_mesh,
    rotation_matrix = full_view,
    cols  = col_AOTU025,
    alpha = 1
  )
p

# ---- Save
ggsave("aotu_morphology.png", plot = p, width = 8, height = 8, dpi = 300)
