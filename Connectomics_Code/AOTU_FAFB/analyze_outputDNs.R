# analyze_outputDNs.R
#
# This script identifies the downstream descending neurons (DNs) for each neuron 
# in a predefined set of cell types and estimates their connection weights. 
# It then analyzes their output projections to VNC neuropils and visualizes the 
# likely contribution of each NOI (Neuron of Interest) to different neuropils.
#
# The script performs the following steps:
#
# 1. **Initialize Parameters**:
#    - Set the working directory and load FlyWire connectivity settings.
#    - Load previously identified `celltypes_of_interest` from a CSV file.
#    - Define directories for plots and data.
#
# 2. **Analyze DN Connectivity for Each Neuron Type**:
#    - Extract connectivity data for right and left hemisphere neurons separately.
#    - Compute synaptic weights for each DN connection.
#    - Calculate the average DN count for each neuron type.
#    - Store DN connectivity results.
#
# Dependencies:
# - Requires `flywire_partner_summary2()`, `flytable_meta()`
#   for connectivity analysis.
# - Uses `ggplot2` for visualization.
#
# Author: MC
# Date: 02/13/2025
#
# Initialize parameters---------------------------------------------------------
# set directory
main_dir = '/Users/mattcollie/Dropbox (HMS)/Connectomics/AOTU' #mac
#main_dir = 'C:/Users/wilson/Dropbox (HMS)/Connectomics/AOTU' #pc
setwd(main_dir)

# fetch settings, data, and libraries
source("initialize_flywire.R")

# set subdirectories
plot_dir <- file.path(main_dir, "plots")
data_dir <- file.path(main_dir, "data")

# Define the input file path in main_dir
input_file <- file.path(data_dir, "celltypes_of_interest.csv")

# Load the saved cell types of interest
if (file.exists(input_file)) {
  celltypes_of_interest <- read.csv(input_file)$cell_type
} else {
  warning("celltypes_of_interest.csv not found. Ensure it has been generated before running this script.")
}

## Plot all steering DNs--------------------------------------------------------
nNOI <- length(celltypes_of_interest)
dn_list <- c('DNa01','DNa02','DNa03','DNae003','DNa11','DNg13','DNg31','DNb02','DNb04','DNb05','DNp09')

# Initialize a list to store results for all cell types
DN_weights <- list()

# Initialize a list to store all DN names encountered
all_DN_names <- list()

# Loop through each cell type
for (i in 1:nNOI) {
  this_celltype <- celltypes_of_interest[i]
  print(paste('Analyzing', this_celltype, '...'))
  
  # Fetch metadata
  this_meta <- flytable_meta(this_celltype)
  this_metaR <- subset(this_meta, side == "right")  # Right side neurons
  this_metaL <- subset(this_meta, side == "left")   # Left side neurons
  
  # Initialize storage for weights and DN counts
  weights_R <- matrix(0, nrow = dim(this_metaR)[1], ncol = length(dn_list))
  weights_L <- matrix(0, nrow = dim(this_metaL)[1], ncol = length(dn_list))
  nDN_R <- numeric(dim(this_metaR)[1])  # Total DN counts for right side neurons
  nDN_L <- numeric(dim(this_metaL)[1])  # Total DN counts for left side neurons
  
  # Store DN names for this cell type
  celltype_DN_names <- c()
  
  # For each neuron on the right side
  for (n in 1:dim(this_metaR)[1]) {
    output_all <- flywire_partner_summary2(this_metaR[n, ], partners = 'out', threshold = min_syn)
    
    # Extract and store DN names
    dn_rows <- output_all[grepl("DN", output_all$cell_type), ]
    if (nrow(dn_rows) > 0) {
      celltype_DN_names <- c(celltype_DN_names, dn_rows$cell_type)
    }
    
    # Count total DNs
    nDN_R[n] <- nrow(dn_rows)
    
    # Store weight values
    for (d in 1:length(dn_list)) {
      thisDN <- dn_list[d]
      output_thisDN <- output_all[grepl(thisDN, output_all$cell_type), ]
      weights_R[n, d] <- sum(output_thisDN$weight, na.rm = TRUE)
    }
  }
  
  # For each neuron on the left side
  for (n in 1:dim(this_metaL)[1]) {
    output_all <- flywire_partner_summary2(this_metaL[n, ], partners = 'out', threshold = min_syn)
    
    # Extract and store DN names
    dn_rows <- output_all[grepl("DN", output_all$cell_type), ]
    if (nrow(dn_rows) > 0) {
      celltype_DN_names <- c(celltype_DN_names, dn_rows$cell_type)
    }
    
    # Count total DNs
    nDN_L[n] <- nrow(dn_rows)
    
    # Store weight values
    for (d in 1:length(dn_list)) {
      thisDN <- dn_list[d]
      output_thisDN <- output_all[grepl(thisDN, output_all$cell_type), ]
      weights_L[n, d] <- sum(output_thisDN$weight, na.rm = TRUE)
    }
  }
  
  # Calculate the mean weights across left and right for each DN
  mean_weights <- colMeans(rbind(weights_R, weights_L), na.rm = TRUE)
  
  # Calculate the average total DN count across left and right, rounded to the nearest whole number
  avg_DN_count <- round(rowMeans(cbind(nDN_R, nDN_L), na.rm = TRUE))
  
  # Store results for this cell type
  DN_weights[[this_celltype]] <- list(
    weights_R = weights_R,
    weights_L = weights_L,
    mean_weights = mean_weights,
    avg_DN_count = avg_DN_count
  )
  
  # Store all DN names found for this cell type
  all_DN_names[[this_celltype]] <- unique(celltype_DN_names)
}

# Define custom cell type and DN orders
desired_order <- c("AOTU019", "AOTU025", "CB0359", "AOTUv3B", "AOTU027",
                   "LAL027", "AOTU015b", "AOTU026", "AOTU015a", "CB2070",
                   "CB3127", "AOTU012", "CB0356", "LAL029")
s
dn_priority <- c("DNa02", "DNa03")  # Desired stacking order: base first
# Re-use the originally defined steering DNs only
dn_list <- c('DNa01','DNa02','DNa03','DNae003','DNa11',
             'DNg13','DNg31','DNb02','DNb04','DNb05','DNp09')

# Stack DNa02 and DNa03 first (base of the bars), others follow
dn_stack_order <- c(setdiff(dn_list, c("DNa02", "DNa03")),"DNa03","DNa02")

# Match actual cell type names to partial patterns and return unique matches
match_order <- function(celltypes, desired_order) {
  matched <- sapply(desired_order, function(pattern) {
    match_idx <- grep(pattern, celltypes, ignore.case = TRUE)
    if (length(match_idx) > 0) {
      return(celltypes[match_idx[1]])
    } else {
      return(NA)
    }
  })
  unique(na.omit(matched))
}

# Build data frame for the stacked bar plot, filling in missing weights with zeros
stacked_weights_df <- data.frame(
  CellType = rep(celltypes_of_interest, each = length(dn_list)),
  DN = rep(dn_list, times = length(celltypes_of_interest)),
  MeanWeight = unlist(lapply(DN_weights, function(x) {
    weights <- rep(0, length(dn_list))
    weights[seq_along(x$mean_weights)] <- x$mean_weights
    return(weights)
  }))
)

# Reorder CellType factor for top-to-bottom plotting
celltype_levels_weights <- rev(match_order(unique(stacked_weights_df$CellType), desired_order))
stacked_weights_df$CellType <- factor(stacked_weights_df$CellType, levels = celltype_levels_weights)

# Set DN stacking order with DNa02 at the base
stacked_weights_df$DN <- factor(stacked_weights_df$DN, levels = dn_stack_order)

# Plot: Stacked horizontal bar of mean weights per DN and cell type
stacked_bar_plot <- ggplot(stacked_weights_df, aes(x = MeanWeight, y = CellType, fill = DN)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Mean Weights for Steering DNs",
    x = "Mean Weight",
    y = "Cell Type"
  ) +
  xlim(0, 500) +  # Set x-axis limit to 500
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1),
    legend.position = "bottom"
  )

# Store the plot
plots <- list()
plots[[1]] <- stacked_bar_plot

# Create DN count data frame
dn_count_df <- data.frame(
  CellType = celltypes_of_interest,
  DNCount = sapply(DN_weights, function(x) mean(x$avg_DN_count, na.rm = TRUE))
)

# Reorder CellType factor similarly for DN count plot
celltype_levels_count <- rev(match_order(unique(dn_count_df$CellType), desired_order))
dn_count_df$CellType <- factor(dn_count_df$CellType, levels = celltype_levels_count)

# Plot: Bar plot of DN counts per cell type
dn_count_plot <- ggplot(dn_count_df, aes(x = DNCount, y = CellType)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    title = "Number of DNs",
    x = "Number of DNs",
    y = "Cell Type"
  ) +
  xlim(0, 500) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(angle = 0, hjust = 1),
    legend.position = "none"
  )

# Add second plot
plots[[2]] <- dn_count_plot

# Save combined figure as SVG and PNG
plotname <- "Stacked_Mean_Weights_and_DNCounts_Horizontal_Corrected"

svglite(paste0(plotname, ".svg"), width = 10, height = 15)
do.call(grid.arrange, c(plots, nrow = 2))
dev.off()

png(paste0(plotname, ".png"), units = "in", width = 10, height = 15, res = 300)
do.call(grid.arrange, c(plots, nrow = 2))
dev.off()
