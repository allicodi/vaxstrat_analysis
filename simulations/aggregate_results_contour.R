
cargs <- commandArgs(trailingOnly = TRUE)
setting <- cargs[[1]]

# Folder with your files
folder <- "~/vegrowth_analysis/simulations/results/contour/"

# List all files in folder (adjust pattern if needed)
files <- list.files(folder, 
		    pattern = paste0("^", setting, "_seed_[0-9]+\\.Rds$"),
		    full.names = TRUE)

# Read each file and store in a list
data_list <- lapply(files, readRDS)

# Combine all data.frames or matrices by row
big_data <- do.call(rbind, data_list)

# Save combined data
saveRDS(big_data, file = file.path(folder, paste0(setting, "_combined_contour_data.Rds")))
