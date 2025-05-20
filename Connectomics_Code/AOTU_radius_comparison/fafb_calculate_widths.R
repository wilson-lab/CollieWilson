# Source packages
source("R/startup/startup.R")
labels <- c(2,3,7,4)

###############################
### generate data on widths ###
###############################

# Functions to estimate neuron cable thickness
estimate_cable_thickness.neuronlist <- function(x, label = NULL, ...){
  nl <- nat::nlapply(x, label=label, FUN = estimate_cable_thickness.neuron, ...)
  nams <- names(nl)
  y <- do.call(rbind, nl)
  y <- as.data.frame(y)
  y$id <- nams
  y
}
estimate_cable_thickness.neuron <- function(x, label =  NULL){
  if(!is.null(label)){
    correct.labels <- x$d$Label == label
    widths <- x$d$W[correct.labels]
  }else{
    label = "all"
    widths <- x$d$W
  }
  widths <- widths[!widths==0]
  med <- median(widths, na.rm = TRUE)
  IQR <- IQR(widths, na.rm = TRUE)
  data.frame(width_median=med,
             width_iqr=IQR,
             width_label=hemibrainr:::standard_compartments(label))
}

# Calculate radii
split.files <-  list.files('/Users/abates/Dropbox (HMS)/hemibrainr/flywire_neurons/630/split', full.names = TRUE)
# list.files("data/swc/", full.names = TRUE)

# For each compartment type....
for(label in labels){

  # Set up parallel backend (adjust ncores as needed)
  doParallel::registerDoParallel(cores = 6)  # You can change the number of cores to use

  # Create progress bar
  pb <- txtProgressBar(max=length(split.files), style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)

  # Parallel for loop using foreach with progress bar
  result_list <- foreach(i = 1:length(split.files),
                         .combine = rbind,
                         .inorder = FALSE,
                         .options.snow=opts) %dopar% {
                           try({
                             x <- read.neurons(split.files[i])
                             thick <- estimate_cable_thickness.neuronlist(x, label = label)
                             thick
                           })
                         }

  # Stop parallel backend
  doParallel::stopImplicitCluster()

  # Save this information elsewhere
  readr::write_csv(result_list, file = sprintf("data/fw_widths_%s.csv",snakecase::to_snake_case(hemibrainr:::standard_compartments(label))))
}
