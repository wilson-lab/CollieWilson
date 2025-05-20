## INITIALIZE ------------------------------------------------------------------
# set working directory (adjust for Mac or PC)
#main_dir = "/Users/mattcollie/Dropbox (HMS)/Connectomics/BANC" # Mac
main_dir = "C:/Users/wilson/HMS Dropbox/Matt Collie/Connectomics/BANC" # PC
setwd(main_dir)

# load necessary libraries
library(bancr)
library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(bit64)
library(ggrepel)
library(writexl)
library(igraph)
source(file.path(main_dir, "getOutputNeurons.R"))

# get BANC edgelist and meta data from 05/06/2025
load("250506_el_meta.RData")

# set subdirectories
plot_dir <- file.path(main_dir, "plots")
# Check if the directory exists, if not, create it
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

# filter banc_meta, eliminate tagged items
filt_banc_meta <- banc_meta %>%
  filter(!grepl("DELETE|MERGE_MONSTER|GLIA|TRACHEA|NOT_A_NEURON",status))

## FETCH LC10a to DN pathway ---------------------------------------------------
# Filter LC10a cells
LC10a_meta <- banc_meta %>%
  filter(str_detect(cell_type, paste('LC10a', collapse = "|"))) %>%
  filter(!is.na(nucleus_supervoxel_id)) %>%
  mutate(
    root_id = as.integer64(root_id),
    cell_type = str_remove(cell_type, "^auto:")  # Remove prefix
  )

# Fetch outputs from each LC10a cell
all_outputs <- LC10a_meta %>%
  select(root_id, cell_type, side) %>%
  mutate(outputs = map(root_id, ~getOutputNeurons(.x, el, filt_banc_meta))) %>%
  unnest(outputs, names_sep = ".") %>%
  distinct(outputs.post_pt_root_id, .keep_all = TRUE)

minWeight <- 10

# Initialize results list
dn_counts_from_LC10a_targets <- data.frame(
  downstream_cell_id = character(),
  downstream_cell_type = character(),
  num_DN_targets = integer(),
  stringsAsFactors = FALSE
)

# Loop through each row in all_outputs
for (i in 1:nrow(all_outputs)) {
  # Get downstream cell ID and cleaned cell type
  post_id <- all_outputs$outputs.post_pt_root_id[i]
  cell_type <- stringr::str_remove(all_outputs$outputs.cell_type[i], "^auto:")
  
  # Run getOutputNeurons on that downstream cell
  outputs <- getOutputNeurons(post_id, el, filt_banc_meta)
  
  # Filter for descending outputs with weight > minWeight
  dn_targets <- outputs %>%
    filter(n > minWeight, cell_class == "descending") %>%
    distinct(post_pt_root_id)
  
  # Store result
  dn_counts_from_LC10a_targets[i, ] <- list(
    downstream_cell_id = post_id,
    downstream_cell_type = cell_type,
    num_DN_targets = nrow(dn_targets)
  )
}
dn_counts_filtered <- dn_counts_from_LC10a_targets %>%
  filter(num_DN_targets > 0) %>%
  drop_na()


## FETCH AOTU OUTPUT DATA ------------------------------------------------------
# Define AOTU cell types of interest
AOTU_cells <- c('AOTU019', 'AOTU025')

# Filter AOTU cells and clean up names
AOTU_meta <- banc_meta %>%
  filter(str_detect(cell_type, paste(AOTU_cells, collapse = "|"))) %>%
  filter(!is.na(nucleus_supervoxel_id)) %>%
  mutate(
    root_id = as.integer64(root_id),
    cell_type = str_remove(cell_type, "^aotu:")  # Remove "aotu:" prefix
  )

# Fetch outputs from each AOTU cell
all_outputs <- AOTU_meta %>%
  select(root_id, cell_type, side) %>%
  mutate(outputs = map(root_id, ~getOutputNeurons(.x, el, filt_banc_meta))) %>%
  unnest(outputs, names_sep = ".")

# Filter for descending neurons
dn_outputs <- all_outputs %>%
  filter(outputs.cell_class == "descending")

# Add AOTU type (AOTU019 vs AOTU025) from presynaptic cell_type
dn_outputs <- dn_outputs %>%
  mutate(AOTU_type = str_extract(cell_type, "AOTU0\\d+"))

# Average output weight per AOTU (left/right) to each DN cell_type
dn_summary <- dn_outputs %>%
  group_by(AOTU_type, side, outputs.cell_type) %>%
  summarise(avg_weight = mean(outputs.n, na.rm = TRUE), .groups = "drop")

# Average across left/right copies for each DN and AOTU cell type
dn_final <- dn_summary %>%
  group_by(outputs.cell_type, AOTU_type) %>%
  summarise(mean_weight = mean(avg_weight, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = AOTU_type, values_from = mean_weight)

dn_final <- dn_final %>%
  rename(
    DN_cell_type = outputs.cell_type,
    AOTU019_weight = AOTU019,
    AOTU025_weight = AOTU025
  ) %>%
  mutate(Mean_weight = rowMeans(select(., AOTU019_weight, AOTU025_weight), na.rm = TRUE)) %>%
  arrange(desc(Mean_weight))

# Check and enforce uniqueness of DN_cell_type in dn_final
dn_final <- dn_final %>%
  distinct(DN_cell_type, .keep_all = TRUE)

# Load and filter functional class table
dn_functional_classes <- read.csv("banc_dn_functional_classes_by_neuron.csv", stringsAsFactors = FALSE)

# Join only the cluster column
dn_final <- dn_final %>%
  left_join(
    dn_functional_classes %>%
      distinct(cell_type, .keep_all = TRUE) %>%
      select(cell_type, cluster),
    by = c("DN_cell_type" = "cell_type")
  )

## DIRECTED PLOT ---------------------------------------------------------------
# Step 1: Filter and reshape edge data
minWeight <- 15

edges <- dn_final %>%
  filter(AOTU019_weight >= minWeight | AOTU025_weight >= minWeight) %>%
  pivot_longer(
    cols = c(AOTU019_weight, AOTU025_weight),
    names_to = "AOTU_source",
    values_to = "weight"
  ) %>%
  filter(weight >= minWeight) %>%
  mutate(
    AOTU_source = ifelse(AOTU_source == "AOTU019_weight", "AOTU019", "AOTU025")
  ) %>%
  select(from = AOTU_source, to = DN_cell_type, weight) %>%
  left_join(dn_final %>% select(DN_cell_type, cluster), by = c("to" = "DN_cell_type"))

# Step 2: Order DNs by cluster
# Set cluster order manually, with DN_20 first
custom_cluster_order <- c(
  "DN_10", "DN_11", "DN_08", "DN_06", "DN_20",
  setdiff(sort(unique(edges$cluster)), c("DN_10", "DN_11", "DN_08", "DN_06", "DN_20"))
)
dn_order <- edges %>%
  distinct(to, cluster) %>%
  mutate(cluster = factor(cluster, levels = custom_cluster_order)) %>%
  arrange(cluster, to) %>%
  pull(to)

# Step 3: Create igraph object
g <- graph_from_data_frame(edges, directed = TRUE)

# Step 4: Assign layout
n_dn <- length(dn_order)
x_dn <- seq(1, by = 5, length.out = n_dn)
x_inputs <- c(mean(range(x_dn)) - 2, mean(range(x_dn)) + 2)

layout <- matrix(nrow = length(V(g)), ncol = 2)
vertex_names <- V(g)$name

layout[vertex_names == "AOTU019", ] <- c(x_inputs[1], 1)
layout[vertex_names == "AOTU025", ] <- c(x_inputs[2], 1)
layout[match(dn_order, vertex_names), ] <- cbind(x_dn, 0)

# Step 5: Assign colors by cluster
dn_clusters <- edges %>%
  distinct(to, cluster) %>%
  arrange(to)

# Get distinct clusters
unique_clusters <- sort(unique(dn_clusters$cluster))
cluster_palette <- brewer.pal(max(3, length(unique_clusters)), "Set3")
cluster_colors <- setNames(cluster_palette[1:length(unique_clusters)], unique_clusters)

# Assign vertex colors
vertex_colors <- case_when(
  vertex_names == "AOTU019" ~ "lightblue",
  vertex_names == "AOTU025" ~ "plum1",
  TRUE ~ cluster_colors[as.character(dn_clusters$cluster[match(vertex_names, dn_clusters$to)])]
)

# Assign edge colors by AOTU source
edge_sources <- ends(g, E(g))[,1]
edge_colors <- ifelse(edge_sources == "AOTU019", "lightblue", "plum1")

# Step 6: Plot PNG
png(file.path(plot_dir, "aotu_weighted_outputs_clustered.png"), width = 1600, height = 900, res = 150)
plot(
  g,
  layout = layout,
  edge.width = E(g)$weight / max(E(g)$weight) * 10,
  edge.arrow.size = 0.3,
  edge.color = edge_colors,
  vertex.color = vertex_colors,
  vertex.label.cex = 0.4,
  vertex.size = 10,
  asp = 0.5
)
legend("bottomleft", legend = unique_clusters, fill = cluster_colors, title = "Cluster", cex = 0.7, bty = "n")
dev.off()

# Step 7: Plot SVG
svg(file.path(plot_dir, "aotu_weighted_outputs_clustered.svg"), width = 14, height = 8)
plot(
  g,
  layout = layout,
  edge.width = E(g)$weight / max(E(g)$weight) * 10,
  edge.arrow.size = 0.3,
  edge.color = edge_colors,
  vertex.color = vertex_colors,
  vertex.label.cex = 0.4,
  vertex.size = 10,
  asp = 0.5
)
legend("bottomleft", legend = unique_clusters, fill = cluster_colors, title = "Cluster", cex = 0.7, bty = "n")
dev.off()

## ANALYZE OUTPUT NEUROPILS ----------------------------------------------------
# Load VNC neuropil annotations from SR
dn_vnc_neuropils <- read.csv("dn_vnc_neuropils.csv", stringsAsFactors = FALSE, colClasses = c("character", "character", "numeric"))
dn_vnc_neuropils$neuron <- banc_latestid(dn_vnc_neuropils$neuron)

# Fetch DN cell types with weight threshold
dn_celltypes <- dn_outputs %>%
  filter(outputs.n >= minWeight) %>%
  select(
    cell_type = outputs.cell_type,
    root_id = outputs.post_pt_root_id
  ) %>%
  distinct(root_id, .keep_all = TRUE)

# Ensure all unique neuropils are added as columns in dn_celltypes
unique_neuropils <- unique(dn_vnc_neuropils$neuropil)
for (neuropil in unique_neuropils) {
  dn_celltypes[[neuropil]] <- NaN
}

# Fill neuropil columns with summed counts
for (i in 1:nrow(dn_celltypes)) {
  this_ID <- dn_celltypes$root_id[i]
  matching_rows <- dn_vnc_neuropils[dn_vnc_neuropils$neuron == this_ID, ]
  
  if (nrow(matching_rows) > 0) {
    for (neuropil in unique_neuropils) {
      total_count <- sum(matching_rows$count[matching_rows$neuropil == neuropil], na.rm = TRUE)
      dn_celltypes[i, neuropil] <- total_count
    }
  }
}

# Remove empty columns (all NA or 0)
empty_cols <- sapply(dn_celltypes, function(col) all(is.na(col) | col == 0))
dn_celltypes <- dn_celltypes[, !empty_cols]

# Save as Excel
write_xlsx(dn_celltypes, path = file.path(plot_dir, "dn_celltypes.xlsx"))

# Summarize across COURT_vnc_ columns
dn_counts <- dn_celltypes %>%
  group_by(cell_type) %>%
  summarise(across(starts_with("COURT_vnc_"), ~ mean(.x, na.rm = TRUE))) %>%
  ungroup()

# Pivot longer to tidy format
dn_long <- dn_counts %>%
  pivot_longer(
    cols = starts_with("COURT_vnc_"),
    names_to = "neuropil",
    values_to = "mean_output"
  ) %>%
  mutate(neuropil = str_remove(neuropil, "^COURT_vnc_"))

# Join cluster info from dn_final
dn_long <- dn_long %>%
  left_join(dn_final %>% select(DN_cell_type, cluster), by = c("cell_type" = "DN_cell_type"))

# Arrange cell types by cluster and fix factor order
dn_long <- dn_long %>%
  mutate(cluster = factor(cluster, levels = custom_cluster_order)) %>%
  arrange(cluster, cell_type) %>%
  mutate(cell_type = factor(cell_type, levels = unique(cell_type)))

# Define priority neuropils and order others
priority_neuropils <- c("ProNM-T1", "MesoNM-T2", "MetaNM-T3")
all_neuropils <- unique(dn_long$neuropil)
priority_present <- intersect(priority_neuropils, all_neuropils)
other_present <- setdiff(all_neuropils, priority_present)

dn_long <- dn_long %>%
  mutate(neuropil = factor(neuropil, levels = rev(c(priority_present, sort(other_present)))))

# Plot: horizontal stacked bars organized by cluster
dn_plot <- ggplot(dn_long, aes(y = cell_type, x = mean_output, fill = neuropil)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ cluster, scales = "free_y", ncol = 1, strip.position = "right") +
  labs(
    title = "Mean Output per Neuropil by DN Cell Type (Clustered)",
    y = "DN Cell Type",
    x = "Mean Output Count",
    fill = "Neuropil"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    axis.text.y = element_text(size = 7),
    legend.position = "right"
  )

# Save plot
ggsave(file.path(plot_dir, "dn_vnc_neuropils_clustered.png"),
       plot = dn_plot, width = 10, height = 12, dpi = 300, bg = "white")

ggsave(file.path(plot_dir, "dn_vnc_neuropils_clustered.svg"),
       plot = dn_plot, width = 10, height = 12, dpi = 300, device = "svg", bg = "white")
