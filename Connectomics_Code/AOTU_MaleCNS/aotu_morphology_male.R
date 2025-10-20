## ───────────────────────────── 0) INITIALIZE ────────────────────────────────

# Optional clean slate
rm(list = ls())

# Paths
main_dir <- "/Users/mattcollie/HMS Dropbox/Matt Collie/Connectomics/AOTU_MaleCNS"
plot_dir <- file.path(main_dir, "plots")
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

# Set colors
col_AOTU019 <- "#75a1e5"
col_AOTU025 <- "#88518f"

# Get region surfaces
brain_mesh <- neuprint_ROI_mesh('CentralBrain',dataset = "male-cns:v0.9")
lo_mesh_l <- neuprint_ROI_mesh('LO(L)',dataset = "male-cns:v0.9")
lo_mesh_r <- neuprint_ROI_mesh('LO(R)',dataset = "male-cns:v0.9")
aotu_mesh_l <- neuprint_ROI_mesh('AOTU(L)',dataset = "male-cns:v0.9")
aotu_mesh_r <- neuprint_ROI_mesh('AOTU(R)',dataset = "male-cns:v0.9")
lal_mesh_l <- neuprint_ROI_mesh('LAL(L)',dataset = "male-cns:v0.9")
lal_mesh_r <- neuprint_ROI_mesh('LAL(R)',dataset = "male-cns:v0.9")

# View matrix (front view)
full_view <- structure(c(0.999475896, -0.002692447, -0.03225674, 0,
                       -0.004669891, -0.998103201, -0.06138369, 0,
                       -0.032030359, 0.061502226, -0.99759269, 0,
                       0,         0,         0,         1),
                     .Dim = c(4L, 4L))

## ───────────────────────── 1) PLOT AOTU019/025  ────────────────────────

# metadata
aotu019_meta <- neuprint_get_meta("AOTU019")
aotu019_meta <- aotu019_meta[grepl("_R", aotu019_meta$name), ]

aotu025_meta <- neuprint_get_meta("AOTU025")
aotu025_meta <- aotu025_meta[grepl("_L", aotu025_meta$name), ]

# Fetch meshes
aotu019_mesh <- neuprint_read_skeletons(aotu019_meta$bodyid,heal.threshold = 20e3, soma = TRUE)
aotu025_mesh <- neuprint_read_skeletons(aotu025_meta$bodyid,heal.threshold = 20e3, soma = TRUE)

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
    root = 0,
    alpha = 1
  )+
  geom_neuron(
    aotu025_mesh,
    rotation_matrix = full_view,
    root = 0,
    cols  = col_AOTU025,
    alpha = 1
  )
p

# ---- Save
setwd(plot_dir)
ggsave("aotu_morphology.png", plot = p, width = 8, height = 8, dpi = 300)


## ───────────────────────── 2) PLOT MINOR AOTU  ────────────────────────
minorAOTU <- c('AOTU002_a','AOTU002_c','AOTU005','AOTU015','AOTU016_a','AOTU016_c',
               'AOTU026','AOTU027','LAL026_a','LAL027','LAL028','LAL029e')
minorColors <- c('#852146','#a72f4b','#c3514c','#d87149','#e49d5e','#edc378',
                 '#f2d49e','#e1e5a1','#b8d285','#9cc155','#648628','#609e8a')

# --- Fetch right-hemisphere skeletons for each minor pathway ---
get_mesh <- function(ct) {
  meta <- tryCatch(neuprint_get_meta(ct), error = function(e) NULL)
  if (is.null(meta) || !"name" %in% names(meta)) return(NULL)
  meta_R <- meta[grepl("_L", meta$name), , drop = FALSE]
  if (nrow(meta_R) == 0) return(NULL)
  tryCatch(
    neuprint_read_skeletons(meta_R$bodyid, heal.threshold = 20e3, soma = TRUE),
    error = function(e) NULL
  )
}

minor_meshes <- setNames(lapply(minorAOTU, get_mesh), minorAOTU)

# Optional: let yourself know if any were missing
missing_types <- names(minor_meshes)[vapply(minor_meshes, is.null, TRUE)]
if (length(missing_types)) message("No _R cells found for: ", paste(missing_types, collapse = ", "))

# --- Plot: same base layers as your original figure ---
roi_col  <- "grey80"
roi_alpha <- 0.3

p_minor <- gganat +
  geom_neuron(
    brain_mesh,
    rotation_matrix = full_view,
    cols  = "grey90",
    alpha = 0.3
  ) +
  geom_neuron(aotu_mesh_l, rotation_matrix = full_view, cols = roi_col, alpha = roi_alpha) +
  geom_neuron(aotu_mesh_r, rotation_matrix = full_view, cols = roi_col, alpha = roi_alpha) +
  geom_neuron(lal_mesh_l,  rotation_matrix = full_view, cols = roi_col, alpha = roi_alpha) +
  geom_neuron(lal_mesh_r,  rotation_matrix = full_view, cols = roi_col, alpha = roi_alpha)

# Add each minor pathway (only those that exist), colored as requested
for (i in seq_along(minor_meshes)) {
  mm <- minor_meshes[[i]]
  if (is.null(mm)) next
  p_minor <- p_minor +
    geom_neuron(
      mm,
      rotation_matrix = full_view,
      root  = 0,
      cols  = minorColors[i],
      alpha = 1
    )
}

# Show and save
p_minor
setwd(plot_dir)
ggsave("aotu_minor_morphology_right.png", plot = p_minor, width = 8, height = 8, dpi = 300)
