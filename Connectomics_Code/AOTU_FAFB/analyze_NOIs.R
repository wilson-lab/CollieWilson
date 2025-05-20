# analyze_NOIs.R
#
# This script identifies neurons downstream of LC10a neurons and estimates their
# contributions to the pursuit circuit by analyzing synaptic connectivity data.
# 
# The script performs the following steps:
# 1. **Initialize parameters**:
#    - Set the main directory and load dependencies.
#    - Create the 'plots' directory if it does not exist.
# 
# 2. **Find downstream cell types of LC10a**:
#    - Fetch metadata for all LC10a neurons.
#    - Identify and store the unique downstream cell types.
#    - Filter out unwanted cell types.
#
# 3. **Compare LC10a input weights for identified cell types**:
#    - Iterate through the selected downstream neurons.
#    - Extract synaptic input from LC10a and output to DNa02 for each neuron.
#    - Compute and store average weights across neurons.
#
# 4. **Filter and visualize results**:
#    - Filter neurons based on minimum synaptic weight thresholds.
#    - Generate bar plots showing LC10a input and DNa02 output weights.
#    - Save the plots as SVG and PNG formats.
#
# 5. **Save identified cell types of interest**:
#    - Export the filtered cell type list as a CSV file in the main directory.
#
# Dependencies:
# - Requires `flywire_partner_summary2()` and `flytable_meta()` functions.
# - Uses `ggplot2`, `gridExtra`, and `svglite` for visualization.
#
# Author: MC 02/12/2025
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
# Check if the directory exists, if not, create it
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}
if (!dir.exists(data_dir)) {
  dir.create(data_dir, recursive = TRUE)
}

# Plot Morphology of Cell Types of Interest-------------------------------------
print('Plotting AOTU neurons...')

# Set working directory for saving plots
setwd(plot_dir)

# Define neuron IDs for AOTU019 and AOTU025
aotu019_id = "720575940631517251"
aotu025_id = "720575940639182424"

# Fetch 3D mesh data for AOTU019 and AOTU025 from FlyWire
aotu019_mesh = read_cloudvolume_meshes(aotu019_id)
aotu025_mesh = read_cloudvolume_meshes(aotu025_id)

# Open a 3D visualization window with a specified view
nopen3d(windowRect = c(2000, 40, 3800, 1000),
        zoom = 0.46,
        userMatrix = structure(c(0.999359846, 0.003123648, -0.035385750, 0,
                                 0.01130394, -0.97230017, -0.23342308, 0,
                                 -0.03367754, -0.23367476, -0.97171968, 0,
                                 0, 0, 0, 1), .Dim = c(4L, 4L)))

# Plot AOTU019 and AOTU025 with their respective colors
plot3d(aotu019_mesh, col = "#81ADB5")  # AOTU019 in blue-green
plot3d(aotu025_mesh, col = "#8D6FA9")  # AOTU025 in purple

# Plot reference structures
plot3d(FAFB14, alpha = vol_o, col = "grey")  # Brain volume
plot3d(FAFB14NP.surf, "LO_L", alpha = glom_o, col = "grey")  # Left Lobula
plot3d(FAFB14NP.surf, "LO_R", alpha = glom_o, col = "grey")  # Right Lobula
plot3d(FAFB14NP.surf, "AOTU_L", alpha = glom_o, col = "grey")  # Left AOTU
plot3d(FAFB14NP.surf, "AOTU_R", alpha = glom_o, col = "grey")  # Right AOTU
plot3d(FAFB14NP.surf, "LAL_L", alpha = glom_o, col = "grey")  # Left LAL
plot3d(FAFB14NP.surf, "LAL_R", alpha = glom_o, col = "grey")  # Right LAL

# Save the visualization as an image
rgl.snapshot("aotu_morphology.png", fmt = 'png')

# Find cell types downstream of LC10a-------------------------------------------
print('Finding cell types downstream of LC10a...')

# fetch all LC10a neurons
lc10_meta = flytable_meta('LC10a*')
nLC = dim(lc10_meta)[1]

# initialize
output_lc10_celltypes = c()

# for each LC10a neuron, find downstream cell types
for (lc in 1:nLC){
  # pull outputs for this LC10a neuron and store cell types
  output_lc10 = flywire_partner_summary2(lc10_meta[lc,], partners = 'out', threshold = min_syn)
  output_lc10_celltypes = c(output_lc10_celltypes,output_lc10$cell_type)
}

# pull unique cell types from output list
output_lc10_celltypes = unique(output_lc10_celltypes)
# omit cell types NOT of interest
cells_to_omit = c('Y3','Y4','LC','1343403608')
print(c('Omitting:',cells_to_omit))
for (o in 1:length(cells_to_omit)){
  output_lc10_celltypes = output_lc10_celltypes[!grepl(cells_to_omit[o],output_lc10_celltypes)]
}
# omit any empty rows
output_lc10_celltypes = output_lc10_celltypes[!is.na(output_lc10_celltypes)]
print('Complete!')

# Compare LC10a input weights for cell types of interest------------------------
print('Comparing LC10a inputs to identified cell types...')
# Set the working directory
setwd(plot_dir)
nNOI = length(output_lc10_celltypes)

# Initialize combined weight containers for LC10a input and DNa02 output weights
lc10_a02_weightsR <- list()
lc10_a02_weightsL <- list()
total_weights <- matrix(NA, nrow = nNOI, ncol = 2)  # to store averages across all neurons for each cell type

# For each cell type of interest
for (i in 1:nNOI) {
  # Pull this cell type
  this_celltype <- output_lc10_celltypes[i]
  print(paste('Analyzing', this_celltype, '...'))
  
  # Fetch metadata
  this_meta <- flytable_meta(this_celltype)
  this_metaR <- subset(this_meta, side == "right")  # right only
  this_metaL <- subset(this_meta, side == "left")   # left only
  
  # Initialize data frames for right and left neurons with two columns: LC10a input and DNa02 output
  lc10_a02_weightsR[[i]] <- data.frame(LC10a_input = numeric(dim(this_metaR)[1]),
                                       DNa02_output = numeric(dim(this_metaR)[1]))
  lc10_a02_weightsL[[i]] <- data.frame(LC10a_input = numeric(dim(this_metaL)[1]),
                                       DNa02_output = numeric(dim(this_metaL)[1]))
  
  # For each neuron of this cell type (right side)
  for (n in 1:dim(this_metaR)[1]) {
    # Fetch all LC10a inputs
    input_all <- flywire_partner_summary2(this_metaR[n, ], partners = 'in', threshold = min_syn)
    input_lc10 <- input_all[grepl("LC10a", input_all$cell_type), ]
    lc10_a02_weightsR[[i]]$LC10a_input[n] <- sum(input_lc10$weight)
    
    # Fetch all DNa02 outputs
    output_all <- flywire_partner_summary2(this_metaR[n, ], partners = 'out', threshold = min_syn)
    output_a02 <- output_all[grepl("DNa02", output_all$cell_type), ]
    lc10_a02_weightsR[[i]]$DNa02_output[n] <- sum(output_a02$weight)
  }
  
  # For each neuron of this cell type (left side)
  for (n in 1:dim(this_metaL)[1]) {
    # Fetch all LC10a inputs
    input_all <- flywire_partner_summary2(this_metaL[n, ], partners = 'in', threshold = min_syn)
    input_lc10 <- input_all[grepl("LC10a", input_all$cell_type), ]
    lc10_a02_weightsL[[i]]$LC10a_input[n] <- sum(input_lc10$weight)
    
    # Fetch all DNa02 outputs
    output_all <- flywire_partner_summary2(this_metaL[n, ], partners = 'out', threshold = min_syn)
    output_a02 <- output_all[grepl("DNa02", output_all$cell_type), ]
    lc10_a02_weightsL[[i]]$DNa02_output[n] <- sum(output_a02$weight)
  }
  
  # Calculate the mean weight across left and right for each cell type
  total_weights[i, 1] <- mean(c(lc10_a02_weightsR[[i]]$LC10a_input, lc10_a02_weightsL[[i]]$LC10a_input), na.rm = TRUE)
  total_weights[i, 2] <- mean(c(lc10_a02_weightsR[[i]]$DNa02_output, lc10_a02_weightsL[[i]]$DNa02_output), na.rm = TRUE)
}

# Prepare data for plotting
minWeight = 10;
valid_indices <- which(total_weights[, 1] > minWeight & total_weights[, 2] > minWeight)  # Filter for weights in both LC10a input and DN output
filtered_celltypes <- output_lc10_celltypes[valid_indices]  # Filter the cell types
filtered_weights <- total_weights[valid_indices, ]  # Filter the weights

# Create filtered data frames
lc10a_data <- data.frame(Neuron = filtered_celltypes, Weight = filtered_weights[, 1])
dna02_data <- data.frame(Neuron = filtered_celltypes, Weight = filtered_weights[, 2])

# Set plot name
plotname <- "LC10a_and_DNa02_Weights"

# Create the plots as in previous code
lc10a_plot <- ggplot(lc10a_data, aes(x = reorder(Neuron, -Weight), y = Weight)) +
  geom_bar(stat = "identity", fill = "firebrick") +
  labs(x = "Neuron Index", y = "LC10a Input Weight", title = "LC10a Input Weights (Sorted by Strength)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dna02_plot <- ggplot(dna02_data, aes(x = reorder(Neuron, -lc10a_data$Weight), y = Weight)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Neuron Index", y = "DNa02 Output Weight", title = "DNa02 Output Weights (Sorted by LC10a Order)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plots as SVG
svglite(paste(plotname, "svg", sep = "."), width = 8, height = 12)
grid.arrange(lc10a_plot, dna02_plot, ncol = 1)
dev.off()

# Save the plots as PNG
png(paste(plotname, "png", sep = "."), units = "in", width = 8, height = 12, res = 1200)
grid.arrange(lc10a_plot, dna02_plot, ncol = 1)
dev.off()

## Store cell types of interest-------------------------------------------------
celltypes_of_interest <- filtered_celltypes

# Define the output file path in main_dir
output_file <- file.path(data_dir, "celltypes_of_interest.csv")

# Save the list as a CSV file
write.csv(data.frame(cell_type = celltypes_of_interest), output_file, row.names = FALSE)
