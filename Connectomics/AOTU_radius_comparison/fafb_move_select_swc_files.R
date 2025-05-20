# Source packages
source("R/startup/startup.R")

# Selected neurons
fw.meta.select <- fw.meta %>%
  dplyr::filter(grepl("^LC|^AOTU",cell_type)|grepl("^LC|^AOTU",hemibrain_type))
fw.meta.select.ids <- unique(fw.meta.select$root_630)

# Get all split files on remote
split.files <- list.files('/Users/abates/Dropbox (HMS)/hemibrainr/flywire_neurons/630/split', full.names = TRUE)
ids <- gsub("\\.swc","",basename(split.files))
split.files.select <- split.files[ids%in%fw.meta.select.ids]

# Copy to local
save.folder <- "data/swc/"
dir.create(save.folder)

# Copy the file with handling potential errors
file.copy(split.files.select, save.folder, overwrite = TRUE)
