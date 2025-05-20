# analyze_outputDNs_BANC
# ------------------------------------------------------------------------------
# Initialize
# ------------------------------------------------------------------------------
# Set working directory (adjust for Mac or PC)
#main_dir = "/Users/mattcollie/Dropbox (HMS)/Connectomics/AOTU" # Mac
main_dir = "C:/Users/wilson/HMS Dropbox/Matt Collie/Connectomics/AOTU" # PC
setwd(main_dir)

# Load necessary libraries
library(bancr)
library(dplyr)

load("250506_el_meta.RData")

# Set search parameters
min_syn = 5  # Minimum number of synapses for filtering connections

# Get BANC meta data
banc.meta = banc_meta

# Get BANC connectivity
banc.el <- el
banc.edgelist.simple <- banc.el %>%
  dplyr::mutate(pre = as.character(pre_pt_root_id),
                post = as.character(post_pt_root_id)) %>%
  dplyr::rename(count = n) %>%
  dplyr::group_by(pre_pt_root_id) %>%
  dplyr::mutate(pre_count = sum(count, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(post_pt_root_id) %>%
  dplyr::mutate(post_count = sum(count, na.rm = TRUE)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(norm = round(count/post_count,6)) %>%
  dplyr::filter(pre_pt_root_id!=post_pt_root_id) %>%
  dplyr::arrange(dplyr::desc(count)) %>%
  dplyr::distinct()

banc.meta.pre <- banc.meta
colnames(banc.meta.pre) <- paste0("pre_",colnames(banc.meta.pre))
banc.meta.post <- banc.meta
colnames(banc.meta.post) <- paste0("post_",colnames(banc.meta.post))

# Add meta data to list
banc.edgelist.simple <- banc.edgelist.simple %>%
  dplyr::left_join(banc.meta.post  %>% 
                     dplyr::select(post_id, post_top_nt, post_cluster,
                                   post_side, post_region, post_super_class, 
                                   post_hemilineage, post_modality, post_nerve, 
                                   post_cell_class, post_cell_sub_class, post_cell_type, post_composite_cell_type) %>%
                     dplyr::distinct(post_id, .keep_all = TRUE),
                   by = c("post"="post_id")) %>%
  dplyr::left_join(banc.meta.pre %>% 
                     dplyr::select(pre_id, pre_top_nt, pre_cluster,
                                   pre_side, pre_region, pre_super_class, 
                                   pre_hemilineage, pre_modality, pre_nerve, 
                                   pre_cell_class, pre_cell_sub_class, pre_cell_type, pre_composite_cell_type) %>%
                     dplyr::distinct(pre_id, .keep_all = TRUE),
                   by = c("pre"="pre_id"))

# ------------------------------------------------------------------------------
# Fetch downstream DNs for circuit of interest
# ------------------------------------------------------------------------------
circuit_of_interest = c('AOTU019','AOTU025')

# Fetch neurons of interest from data table
bm_noi = bm %>% filter(str_detect(cell_type, str_c(circuit_of_interest, collapse = "|")))

# Initialize empty dataframe to store all DN outputs
dn_outputs_all_runs = tibble()

# For each neuron
for (i in 1:nrow(bm_noi)) {
  # Fetch this neuron
  this_neuron = bm_noi[i,]
  
  # Fetch outputs for this cell
  outputs_all = banc_partner_summary(banc_latestid(this_neuron), partners = "output", threshold = min_syn)

  # Fetch meta data for all outputs
  outputs_all_meta = bm %>% filter(root_id %in% outputs_all$post_id)
  
  # Pull metadata
  outputs_dns = outputs_all_meta %>% filter(str_detect(cell_type, "DN"))
  
  # Store results
  dn_outputs_all_runs = bind_rows(dn_outputs_all_runs, outputs_dns)
}

# Filter out repeats
dn_outputs_all_runs = dn_outputs_all_runs %>% distinct(root_id, .keep_all = TRUE)