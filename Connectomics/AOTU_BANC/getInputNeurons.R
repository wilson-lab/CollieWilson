# Function that returns table of inputs for neurons specified by list of 
#  root_ids

getInputNeurons <- function(targetIDs, el, banc_meta) {
  target_inputs <- el %>%
    filter((post_pt_root_id %in% targetIDs))
  target_inputs_info <- target_inputs %>%
  left_join(
    banc_meta %>% 
#      mutate(root_id = as.integer64(root_id)) %>%
      select(root_id, super_class, cell_class, cell_sub_class, cell_type),
    by = c("pre_pt_root_id" = "root_id")
  )
  return(target_inputs_info)
}