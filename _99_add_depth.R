# Install and load necessary packages
# install.packages("dplyr")
# install.packages("geosphere")
# an AI wonder...

library(tidyverse)
library(geosphere) # For calculating geographic distances
setwd("~/GitHub/krill-phyto")
main_dir = getwd()
data_dir = paste0(main_dir,'/data/')

# --- Configuration --- ########################################################
# Define the paths to your input files
# The  primary data file

data_main <- data.frame(read.csv( paste0(data_dir,"PB_NASC_sizeclass2.csv"), header=TRUE)) 
output_file_path <- paste0(data_dir,"PB_NASC_sizeclass2_with_depth.csv")

# the file with depth
data_depth <- data.frame(read.delim("Tolimieri_et_al_2015_Ecosphere_2km_grid.txt", header=TRUE))

# Define column names for latitude and longitude in both files
# Ensure these match your actual column names
lat_col_file1 <- "LAT"
lon_col_file1 <- "LON"

lat_col_file2 <- "LAT"
lon_col_file2 <- "LON"
depth_col_file2 <- "SRTM_M" # The column in file2 that contains depth

# --- Data Preparation ---
# Ensure latitude and longitude columns are numeric
data_main[[lat_col_file1]] <- as.numeric(data_main[[lat_col_file1]])
data_main[[lon_col_file1]] <- as.numeric(data_main[[lon_col_file1]])

data_depth[[lat_col_file2]] <- as.numeric(data_depth[[lat_col_file2]])
data_depth[[lon_col_file2]] <- as.numeric(data_depth[[lon_col_file2]])
data_depth[[depth_col_file2]] <- as.numeric(data_depth[[depth_col_file2]])

# Remove rows with NA in critical columns after conversion
data_main <- data_main %>%
  filter(!is.na(!!sym(lat_col_file1)) & !is.na(!!sym(lon_col_file1)))
data_depth <- data_depth %>%
  filter(!is.na(!!sym(lat_col_file2)) & !is.na(!!sym(lon_col_file2)) & !is.na(!!sym(depth_col_file2)))

message(paste("File 1 has", nrow(data_main), "valid rows."))
message(paste("File 2 has", nrow(data_depth), "valid rows."))

# --- Find Closest Match and Merge Depth ---
message("Finding closest matches and merging depth...")

# Initialize a vector to store the matched depths
matched_depths <- numeric(nrow(data_main))

# Loop through each row in the main data frame
for (i in 1:nrow(data_main)) {
  # Get latitude and longitude for the current row in data_main
  lat1 <- data_main[[lat_col_file1]][i]
  lon1 <- data_main[[lon_col_file1]][i]
  
  # Create a matrix of coordinates for the current point in data_main
  point1 <- c(lon1, lat1)
  
  # Create a matrix of coordinates for all points in data_depth
  points2 <- as.matrix(data_depth[, c(lon_col_file2, lat_col_file2)])
  
  # Calculate the Haversine distance between point1 and all points in data_depth
  # The distHaversine function expects longitude, then latitude
  distances <- distHaversine(point1, points2) # Distance in meters
  
  # Find the index of the minimum distance
  closest_index <- which.min(distances)
  
  # Get the depth from the closest matching row in data_depth
  matched_depths[i] <- data_depth[[depth_col_file2]][closest_index]
  
  # Optional: Print progress for large files
  if (i %% 100 == 0) {
    message(paste("Processed", i, "of", nrow(data_main), "rows."))
  }
}

# Add the matched depths as a new column to the main data frame
data_main$matched_depth <- matched_depths

message("Depth merging complete.")

# --- Write the result to a new CSV file ---
message(paste("Writing results to:", output_file_path))
write.csv(data_main, output_file_path, row.names = FALSE)

message("Process finished successfully!")
message(paste("Output saved to:", output_file_path))
# Display the first few rows of the modified data frame
print(head(data_main))
