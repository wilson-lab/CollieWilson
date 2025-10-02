# Initialize parameters---------------------------------------------------------
# set directory
#main_dir = '/Users/mattcollie/HMS Dropbox/Matt Collie/Connectomics/AOTU_FAFB' #mac
main_dir = 'C:/Users/wilson/HMS Dropbox/Matt Collie/Connectomics/AOTU_FAFB' #pc
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

# Fetch meta data for aIPg cells-------------------------------------------------
library(bit64)

# aIPg is a secondary celltype label, need to fetch from Codex search results
search_results <- read.csv(
  file.path(data_dir, "search_results_aipg.csv"),
  colClasses = c(root_id = "integer64")
)

# make separate data.frames, fetch meta data
aIPga_df <- flytable_meta(subset(search_results, grepl("aIPga", cell_type))["root_id"])
aIPgb_df <- flytable_meta(subset(search_results, grepl("aIPgb", cell_type))["root_id"])
aIPgc_df <- flytable_meta(subset(search_results, grepl("aIPgc", cell_type))["root_id"])

# Fetch output connectivity-----------------------------------------------------
outputs_a <- flywire_partner_summary2(aIPga_df, partners = "out", threshold = min_syn)
outputs_b <- flywire_partner_summary2(aIPgb_df, partners = "out", threshold = min_syn)
outputs_c <- flywire_partner_summary2(aIPgc_df, partners = "out", threshold = min_syn)

# Average weight per cell_type
avg_outputs_a <- aggregate(weight ~ cell_type, data = outputs_a, FUN = function(x) mean(x, na.rm = TRUE))
avg_outputs_b <- aggregate(weight ~ cell_type, data = outputs_b, FUN = function(x) mean(x, na.rm = TRUE))
avg_outputs_c <- aggregate(weight ~ cell_type, data = outputs_c, FUN = function(x) mean(x, na.rm = TRUE))

# (Optional) sort by descending average weight
avg_outputs_a <- avg_outputs_a[order(-avg_outputs_a$weight), ]
avg_outputs_b <- avg_outputs_b[order(-avg_outputs_b$weight), ]
avg_outputs_c <- avg_outputs_c[order(-avg_outputs_c$weight), ]

# Combine outputs from all subtypes
outputs_all <- rbind(outputs_a, outputs_b, outputs_c)

# Average weight per cell_type across a, b, c
avg_outputs_all <- aggregate(weight ~ cell_type, data = outputs_all, FUN = function(x) mean(x, na.rm = TRUE))

# (Optional) sort descending
avg_outputs_all <- avg_outputs_all[order(-avg_outputs_all$weight), ]


# Find input regions for AOTU019------------------------------------------------
AOTU019_meta = flytable_meta('AOTU019*')

# Fetch and process inputs
right_inputs <- flywire_partner_summary2(AOTU019_meta[1, ], partners = 'in', threshold = min_syn) %>%
  mutate(neuropil_std = clean_neuropil(top_np)) %>%
  group_by(neuropil_std) %>%
  summarise(weight_right = sum(weight), .groups = "drop")

left_inputs <- flywire_partner_summary2(AOTU019_meta[2, ], partners = 'in', threshold = min_syn) %>%
  mutate(neuropil_std = clean_neuropil(top_np)) %>%
  group_by(neuropil_std) %>%
  summarise(weight_left = sum(weight), .groups = "drop")

# Merge and average
combined <- full_join(right_inputs, left_inputs, by = "neuropil_std") %>%
  mutate(
    weight_right = replace_na(weight_right, 0),
    weight_left = replace_na(weight_left, 0),
    weight_avg = (weight_right + weight_left) / 2
  ) %>%
  filter(weight_avg >= 50)

# Create the plot and assign it to a variable
p <- ggplot(combined, aes(x = "", y = weight_avg, fill = neuropil_std)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  theme_void() +
  labs(title = "Average Input Synapse Counts by Neuropil (≥ 50, AOTU019 L & R)", fill = "Neuropil") +
  theme(legend.position = "right")

# Save as SVG
ggsave("AOTU019_npinputs_pie.svg", plot = p, device = svglite::svglite, width = 6, height = 6)
# Save as PNG
ggsave("AOTU019_npinputs_pie.png", plot = p, dpi = 300, width = 6, height = 6)

AOTU025_meta = flytable_meta('AOTU025*')

# Fetch and process inputs
right_inputs <- flywire_partner_summary2(AOTU025_meta[1, ], partners = 'in', threshold = min_syn) %>%
  mutate(neuropil_std = clean_neuropil(top_np)) %>%
  group_by(neuropil_std) %>%
  summarise(weight_right = sum(weight), .groups = "drop")

left_inputs <- flywire_partner_summary2(AOTU025_meta[2, ], partners = 'in', threshold = min_syn) %>%
  mutate(neuropil_std = clean_neuropil(top_np)) %>%
  group_by(neuropil_std) %>%
  summarise(weight_left = sum(weight), .groups = "drop")

# Merge and average
combined <- full_join(right_inputs, left_inputs, by = "neuropil_std") %>%
  mutate(
    weight_right = replace_na(weight_right, 0),
    weight_left = replace_na(weight_left, 0),
    weight_avg = (weight_right + weight_left) / 2
  ) %>%
  filter(weight_avg >= 50)

# Create the plot and assign it to a variable
p <- ggplot(combined, aes(x = "", y = weight_avg, fill = neuropil_std)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  theme_void() +
  labs(title = "Average Input Synapse Counts by Neuropil (≥ 50, AOTU025 L & R)", fill = "Neuropil") +
  theme(legend.position = "right")

# Save as SVG
ggsave("AOTU025_npinputs_pie.svg", plot = p, device = svglite::svglite, width = 6, height = 6)
# Save as PNG
ggsave("AOTU025_npinputs_pie.png", plot = p, dpi = 300, width = 6, height = 6)


## Analyze shared inputs to AOTU019/025-----------------------------------------
min_synshared = 20
min_synalone = 200
# Load input data
right_019 <- flywire_partner_summary2(AOTU019_meta[1, ], partners = 'in', threshold = min_syn)
left_019 <- flywire_partner_summary2(AOTU019_meta[2, ], partners = 'in', threshold = min_syn)
right_025 <- flywire_partner_summary2(AOTU025_meta[1, ], partners = 'in', threshold = min_syn)
left_025 <- flywire_partner_summary2(AOTU025_meta[2, ], partners = 'in', threshold = min_syn)

# summarize total input weights per cell type for each neuron
sum_inputs <- function(df) {
  df %>%
    group_by(cell_type) %>%
    summarize(total_weight = sum(weight), .groups = 'drop')
}

sum_right_019 <- sum_inputs(right_019)
sum_left_019  <- sum_inputs(left_019)
sum_right_025 <- sum_inputs(right_025)
sum_left_025  <- sum_inputs(left_025)

# Add suffixes to distinguish left/right
colnames(sum_right_019)[2] <- "weight_right_019"
colnames(sum_left_019)[2]  <- "weight_left_019"
colnames(sum_right_025)[2] <- "weight_right_025"
colnames(sum_left_025)[2]  <- "weight_left_025"

# Combine left and right for each neuron
avg_019 <- full_join(sum_right_019, sum_left_019, by = "cell_type") %>%
  mutate(mean_019 = rowMeans(select(., starts_with("weight_")), na.rm = TRUE))

avg_025 <- full_join(sum_right_025, sum_left_025, by = "cell_type") %>%
  mutate(mean_025 = rowMeans(select(., starts_with("weight_")), na.rm = TRUE))

# Filter out cell types that contain "LC10"
avg_019 <- avg_019 %>% filter(!grepl("LC10", cell_type))
avg_025 <- avg_025 %>% filter(!grepl("LC10", cell_type))

shared_inputs <- full_join(avg_019[, c("cell_type", "mean_019")],
                           avg_025[, c("cell_type", "mean_025")],
                           by = "cell_type")

# Sort by combined input strength
shared_inputs <- shared_inputs %>%
  mutate(combined = rowMeans(cbind(mean_019, mean_025))) %>%
  arrange(desc(combined))

# View top shared inputs
head(shared_inputs, 10)

# Step 1: Filter shared inputs with >20 to both
shared_inputs <- shared_inputs %>% filter(!is.na(cell_type))
shared_inputs_filtered <- shared_inputs %>%
  filter(mean_019 > min_synshared, mean_025 > min_synshared)

# Step 2: Filter exclusive strong inputs (>50 to either) not in shared list
exclusive_019 <- shared_inputs %>%
  filter(mean_019 > min_synalone, !(cell_type %in% shared_inputs_filtered$cell_type))

exclusive_025 <- shared_inputs %>%
  filter(mean_025 > min_synalone, !(cell_type %in% shared_inputs_filtered$cell_type))

# Step 3: Create edge list
edges <- data.frame(
  from = c(rep(shared_inputs_filtered$cell_type, 2),
           exclusive_019$cell_type,
           exclusive_025$cell_type),
  to = c(rep(c("AOTU019", "AOTU025"), each = nrow(shared_inputs_filtered)),
         rep("AOTU019", nrow(exclusive_019)),
         rep("AOTU025", nrow(exclusive_025))),
  weight = c(shared_inputs_filtered$mean_019,
             shared_inputs_filtered$mean_025,
             exclusive_019$mean_019,
             exclusive_025$mean_025)
)

# Step 4: Build graph
g <- graph_from_data_frame(edges, directed = TRUE)

# Step 5: Define node groups
input_shared <- shared_inputs_filtered$cell_type
# Sort AOTU019-only inputs by strength (descending)
exclusive_019 <- exclusive_019 %>%
  arrange(desc(mean_019))
input_019 <- exclusive_019$cell_type
input_025 <- exclusive_025$cell_type
output_nodes <- c("AOTU019", "AOTU025")

# Combine all nodes and ensure uniqueness
all_nodes <- unique(c(input_shared, input_019, input_025, output_nodes))

# Step 6: Create layout
layout <- matrix(NA, nrow = length(all_nodes), ncol = 2)

# Shared inputs (top row)
layout[1:length(input_shared), 1] <- seq(-5, 5, length.out = length(input_shared))
layout[1:length(input_shared), 2] <- 1

# Exclusive AOTU019 (bottom left)
start_019 <- length(input_shared) + 1
end_019 <- start_019 + length(input_019) - 1
layout[start_019:end_019, 1] <- seq(-6, -1, length.out = length(input_019))
layout[start_019:end_019, 2] <- -1

# Exclusive AOTU025 (bottom right) only if not empty
if (length(input_025) > 0) {
  start_025 <- end_019 + 1
  end_025 <- start_025 + length(input_025) - 1
  layout[start_025:end_025, 1] <- seq(1, 3, length.out = length(input_025))
  layout[start_025:end_025, 2] <- -1
} else {
  end_025 <- end_019  # continue from previous end
}

# Output nodes (middle row)
layout[(end_025 + 1), ] <- c(-1, 0)  # AOTU019
layout[(end_025 + 2), ] <- c(1, 0)   # AOTU025

# Safety check: ensure all V(g)$name are in all_nodes
missing_nodes <- setdiff(V(g)$name, all_nodes)
if (length(missing_nodes) > 0) {
  stop("Missing node(s) in layout: ", paste(missing_nodes, collapse = ", "))
}

# Reorder layout to match graph vertex order
layout_order <- match(V(g)$name, all_nodes)
if (any(is.na(layout_order))) stop("Layout contains NA indices — check node matching.")
layout <- layout[layout_order, ]

# Step 7: Plot
# Set output file names
png_file <- file.path(plot_dir, "aotu_inputs_plot.png")
svg_file <- file.path(plot_dir, "aotu_inputs_plot.svg")

# Save as PNG
png(png_file, width = 1200, height = 800, res = 150)
plot(
  g,
  layout = layout,
  edge.width = E(g)$weight / max(E(g)$weight) * 10,
  edge.arrow.size = 0.3,
  vertex.color = case_when(
    V(g)$name %in% output_nodes ~ "orange",
    V(g)$name %in% input_019 ~ "lightgreen",
    V(g)$name %in% input_025 ~ "lightgreen",
    TRUE ~ "skyblue"
  ),
  vertex.label.cex = 0.8,
  vertex.size = 20,
  asp = 0.5
)
dev.off()

# Save as SVG
svg(svg_file, width = 10, height = 7)
plot(
  g,
  layout = layout,
  edge.width = E(g)$weight / max(E(g)$weight) * 10,
  edge.arrow.size = 0.3,
  vertex.color = case_when(
    V(g)$name %in% output_nodes ~ "orange",
    V(g)$name %in% input_019 ~ "lightgreen",
    V(g)$name %in% input_025 ~ "lightgreen",
    TRUE ~ "skyblue"
  ),
  vertex.label.cex = 0.8,
  vertex.size = 20,
  asp = 0.5
)
dev.off()