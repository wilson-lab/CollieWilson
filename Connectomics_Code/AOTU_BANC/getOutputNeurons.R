# Function that returns table of outputs for neurons specified by list of 
#  root_ids

getOutputNeurons <- function(targetIDs, el, banc_meta) {
  target_outputs <- el %>%
    filter(pre_pt_root_id %in% targetIDs)
  
  target_outputs_info <- target_outputs %>%
    left_join(
      banc_meta %>%
        mutate(root_id = as.integer64(root_id)) %>%  # force correct type
        select(root_id, super_class, cell_class, cell_sub_class, cell_type),
      by = c("post_pt_root_id" = "root_id")
    )
  
  return(target_outputs_info)
}