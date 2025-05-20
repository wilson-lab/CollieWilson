# initialize_flywire.R
#
# This script initializes the required packages, sets up access to FlyWire data, 
# and defines essential parameters and visualization settings for downstream analyses.
#
# The script performs the following steps:
#
# 1. **Initialize Packages**:
#    - Loads necessary R libraries (`fafbseg`, `elmr`, `dplyr`, etc.).
#    - Ensures `natmanager` is available to install missing dependencies.
#
# 2. **Fetch FlyWire Data**:
#    - Sets the FlyWire access token for authentication.
#    - Downloads FlyWire release data, including connectivity and cell-type information.
#    - Defines minimum synapse thresholds for different neuron types (LC10, DN).
#
# 3. **Set Plot Variables**:
#    - Defines a pastel color palette for visualizations.
#    - Creates a function (`map2color`) for mapping values to colors.
#    - Generates color scales for heatmaps.
#    - Sets transparency levels for neuropil and glomerulus volumes.
#
# Author: MC 02/12/2025
#
# Initialize Packages-----------------------------------------------------------
# install natmanager if required
#if (!requireNamespace("natmanager")) install.packages("natmanager")
#natmanager::install(pkgs="fafbseg")

# load libraries
library(fafbseg)
library(elmr)
library(dplyr)
library(svglite)
library(RColorBrewer)
library(writexl)
library(ggplot2)
library(gridExtra)
library(reshape2)

# Fetch Flywire Data------------------------------------------------------------
# Record access tokens for FlyWire and Banc
#flywire_set_token()

# Fetch FlyWire release data (connectivity and cell type information)
download_flywire_release_data('all')

# min synapse threshold to be included for all analyses
min_syn = 5
min_synLC10 = 25
min_synDN = 50

# Set Plot Variables------------------------------------------------------------
# Define a pastel color palette with 16+ distinct but harmonious colors
base_colors <- c("#FFB6C1", "#FFD700", "#87CEFA", "#A2CD5A", "#98FB98", "#DDA0DD", "#F4A460", 
                 "#FA8072", "#00CED1", "#FF69B4", "#CD5C5C", "#4682B4", "#9ACD32",
                 "#E6A8D7", "#FFA07A", "#B0C4DE")

# generate heatmap colors
map2color <- function(x, pal, limits = range(x)){
  pal[findInterval(x, seq(limits[1], limits[2], length.out = length(pal) + 1), 
                   all.inside=TRUE)]
}
myPalette = colorRampPalette(c("white","black"))
mapPal_bw = myPalette(100)
mapPal_heat = heat.colors(100)
vol_o = 0.05 #neuropil volume opacity
glom_o = 0.15 #glomerulus volume opacity
