# analyze_inputLC10s.R
#
# This script analyzes LC10a synaptic inputs onto neurons of interest 
# and generates heatmaps to visualize synaptic strength using 3D 
# reconstructions from FlyWire.
# 
# The script performs the following steps:
#
# 1. **Initialize parameters**:
#    - Set the main directory and load dependencies.
#    - Define and create necessary subdirectories (`plots`, `heatmaps`, `data`).
#    - Load the list of cell types of interest from a CSV file.
#
# 2. **Generate LC10a input heatmaps**:
#    - Loop through each neuron of interest and fetch its metadata.
#    - Identify presynaptic inputs from LC10a neurons in the right hemisphere.
#    - Generate 3D mesh visualizations of LC10a inputs.
#    - Create heatmaps of synaptic input strength in both greyscale and color.
#    - Save the resulting heatmaps as PNG images.
#
# 3. **Plot all LC10a neurons for comparison**:
#    - Fetch 3D reconstructions of all LC10a neurons.
#    - Generate a reference plot of all LC10a neurons in black.
#    - Save the plot as a PNG image for comparison.
#
# Dependencies:
# - Requires `flywire_partner_summary2()`, `flytable_meta()`, and 
#   `read_cloudvolume_meshes()` functions.
# - Uses `rgl` for 3D visualization.
#
# Author: MC 02/14/2025
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
heatmap_dir <- file.path(main_dir, "heatmap")
data_dir <- file.path(main_dir, "data")

# Check if the directory exists, if not, create it
if (!dir.exists(heatmap_dir)) {
  dir.create(heatmap_dir, recursive = TRUE)
}

# Define the input file path in main_dir
input_file <- file.path(data_dir, "celltypes_of_interest.csv")

# Load the saved cell types of interest
if (file.exists(input_file)) {
  celltypes_of_interest <- read.csv(input_file)$cell_type
} else {
  warning("celltypes_of_interest.csv not found. Ensure it has been generated before running this script.")
}

# ------------------------------------------------------------------------------
# Generate LC10a Input Heatmaps by Synapse Weight
# ------------------------------------------------------------------------------
print('Generating LC10a synapse heatmaps... (this may take a while)')

# Set working directory for saving heatmaps
setwd(heatmap_dir)

# Get the number of neurons of interest (NOI)
nNOI = length(celltypes_of_interest)

# Loop through each neuron of interest
for (i in 9:nNOI) {
  # Pull this cell type
  this_celltype <- celltypes_of_interest[i]
  print(paste('Analyzing', this_celltype, '...'))
  
  # Fetch metadata
  this_meta <- flytable_meta(this_celltype)
  this_metaR <- subset(this_meta, side == "right")  # Select only right hemisphere neurons
  
  # Check how many rows exist for this cell type in the right hemisphere
  n_metaR = nrow(this_metaR)
  
  # Iterate over each unique neuron of this cell type
  for (j in 1:n_metaR) {
    this_neuron <- this_metaR[j,]  # Select current neuron instance
    
    print(paste("Processing neuron", this_neuron$root_id, "of", this_celltype))
    
    # Fetch all presynaptic input partners for this neuron
    input_all = flywire_partner_summary2(this_neuron, partners = 'in', threshold = min_syn)
    
    # Filter to include only LC10a neurons on the right hemisphere
    input_lc10 = input_all[grepl("LC10a", input_all$cell_type), ]
    input_lc10 = input_lc10[grepl("right", input_lc10$side), ]
    
    # Count the number of LC10a inputs
    nLC10 = nrow(input_lc10)
    
    # Process neurons only if there are more than 10 LC10a inputs
    if (nLC10 > 10) {
      
      # Check for duplicate input IDs
      duplicate_counts <- table(input_lc10$pre_pt_root_id)
      duplicates_only <- duplicate_counts[duplicate_counts > 1]
      
      if (length(duplicates_only) > 0) {
        print("Warning: Duplicate pre-synaptic root IDs detected.")
        print(duplicates_only)
      }
      
      # Fetch 3D meshes for LC10a input neurons
      lc10_mesh = read_cloudvolume_meshes(input_lc10$pre_pt_root_id, 
                                          cloudvolume.url = "graphene://https://prod.flywire-daf.com/segmentation/table/fly_v31")
      
      # Open a 3D visualization window for greyscale heatmap
      nopen3d(windowRect = c(2000, 40, 2900, 1000),
              zoom = 0.71,
              userMatrix = structure(c(-0.8142193, 0.4023438, 0.4184932, 0,
                                       -0.1945970, -0.8683136, 0.4562103, 0,
                                       0.5469477, 0.2900211, 0.7853076, 0,
                                       0, 0, 0, 1), .Dim = c(4L, 4L)))
      
      # Plot brain structure and LC10a inputs in greyscale
      plot3d(FAFB14NP.surf, "LO_L", alpha = vol_o, col = "grey")
      plot3d(lc10_mesh, col = map2color(input_lc10$weight, mapPal_bw))  # Greyscale heatmap
      
      # Save greyscale heatmap
      plotname = paste(this_neuron$cell_type, this_neuron$root_id, "lc10a", sep = "_")
      rgl.snapshot(paste(plotname, 'grey.png', sep = '_'), fmt = 'png')
      close3d()
      
      # Open a new 3D visualization window for color heatmap
      nopen3d(windowRect = c(2000, 40, 2900, 1000),
              zoom = 0.71,
              userMatrix = structure(c(-0.8142193, 0.4023438, 0.4184932, 0,
                                       -0.1945970, -0.8683136, 0.4562103, 0,
                                       0.5469477, 0.2900211, 0.7853076, 0,
                                       0, 0, 0, 1), .Dim = c(4L, 4L)))
      
      # Plot brain structure and LC10a inputs in color scale
      plot3d(FAFB14NP.surf, "LO_L", alpha = vol_o, col = "grey")
      plot3d(lc10_mesh, col = map2color(input_lc10$weight, mapPal_heat))  # Color heatmap
      
      # Save color heatmap
      rgl.snapshot(paste(plotname, 'color.png', sep = '_'), fmt = 'png')
      close3d()
    }
  }
}


print("LC10a synapse heatmap generation complete!")

# ------------------------------------------------------------------------------
# Plot All LC10a Neurons for Comparison
# ------------------------------------------------------------------------------

print('Plotting ALL LC10a neurons for comparison (this may take a while)...')

# Fetch metadata
this_meta <- flytable_meta('LC10a')
this_metaR <- subset(this_meta, side == "right")  # Select only right hemisphere neurons

# Set working directory for saving plots
setwd(heatmap_dir)

# Fetch 3D meshes for all LC10a neurons
lc10_mesh = read_cloudvolume_meshes(this_metaR$root_id, 
                                    cloudvolume.url = "graphene://https://prod.flywire-daf.com/segmentation/table/fly_v31")

# Open a 3D visualization window for plotting all LC10a neurons
nopen3d(windowRect = c(2000, 40, 2900, 1000),
        zoom = 0.71,
        userMatrix = structure(c(-0.8142193, 0.4023438, 0.4184932, 0,
                                 -0.1945970, -0.8683136, 0.4562103, 0,
                                 0.5469477, 0.2900211, 0.7853076, 0,
                                 0, 0, 0, 1), .Dim = c(4L, 4L)))

# Plot brain structure and all LC10a neurons in black
plot3d(FAFB14NP.surf, "LO_L", alpha = vol_o, col = "grey")
plot3d(lc10_mesh, col = "black")  # All LC10a neurons

# Save plot
plotname = "all_lc10a.png"
rgl.snapshot(plotname, fmt = 'png')

print("All LC10a neuron plotting complete!")
