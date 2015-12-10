# 
# This function converts raw data from the output of ca_snapsSS et al. to one 
#   that the spatialwarnings can understand. 
# 

# Base paths
result_folder <- '/home/alex/work/2014-2015/SpatialStress/ifcam_results/'

convert_dat <- function(upper_branch, lower_branch, newfile) { 
  
  # Convert all matrices to binary_matrices (for loop = less memory used)
  for (i in seq.int(upper_branch[[2]])) { 
    print(i)
    upper_branch[[2]][[i]] <- lapply(upper_branch[[2]][[i]], as.binary_matrix)
    lower_branch[[2]][[i]] <- lapply(lower_branch[[2]][[i]], as.binary_matrix)
  }
  
  save(upper_branch, lower_branch, file = newfile, compress = 'bzip2')
  
  rm(lower_branch, upper_branch)
}

# Convert forestgap 
load(paste0(result_folder,'result_forestgap.rda'), verbose = TRUE)
convert_dat(result_forestgap_upper, result_forestgap_lower, 
            newfile = paste0(result_folder, "result_forestgap_processed.rda" ))
rm(result_forestgap_upper, result_forestgap_lower)
gc()

# Convert predprey 
load(paste0(result_folder,'result_predprey.rda'), verbose = TRUE)
convert_dat(result_predprey_upper, result_predprey_lower, 
            newfile = paste0(result_folder, "result_predprey_processed.rda" ))
rm(result_predprey_upper, result_predprey_lower)
gc()

# Convert musselbed
load(paste0(result_folder,'result_musselbed.rda'), verbose = TRUE)
convert_dat(result_musselbed_upper, result_musselbed_lower, 
            newfile = paste0(result_folder, "result_musselbed_processed.rda" ))
rm(result_musselbed_upper, result_musselbed_lower)
gc()

# Convert grazing
load(paste0(result_folder,'result_grazing.rda'), verbose = TRUE)
convert_dat(result_grazing_upper, result_grazing_lower, 
            newfile = paste0(result_folder, "result_grazing_processed.rda" ))
rm(result_grazing_upper, result_grazing_lower)
gc()

