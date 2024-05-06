#### EXAMPLE LIDR PROCESSING WORKFLOW ####

install.packages("lidR")
install.packages("future")
install.packages("RCSF")
install.packages("sp")
install.packages("raster")

library(lidR)
library(future) # for parallel processing

## create a folders to save outputs 
dir = paste0("G:/MS1/lidR_sample")

dir.create(paste0(dir,"/output"))
dir.create(paste0(dir,"/output/GROUND")) # ground classification
dir.create(paste0(dir,"/output/NORM")) # normalized point cloud
dir.create(paste0(dir,"/output/RASTER")) # rasterized DEM & CHM

## read in the LAS Catalog 
LCT = readLAScatalog(folder = paste0(dir, "/D17a_laz"))

#### GROUND CLASSIFICATION ####
# classify ground/canopy from the filtered point cloud
# using cloth simulation function

## options for reading in the LAS Catalog
opt_chunk_buffer(LCT) = 10 # buffer each tile 10m
opt_laz_compression(LCT) = FALSE # write output as .las 
# Point cloud thinning to speed things up by taking a random point within each 0.1m voxel
# consider using higher/lower resolution voxel depending on point cloud density
opt_filter(LCT) = "-thin_with_voxel 0.1" 
# output naming convention will keep the tiling system
opt_output_files(LCT) = paste0(dir, "\\output\\GROUND\\{*}_GROUND") 
opt_progress(LCT) = TRUE

# set up and execute task in parallel (2 cores) 
plan(multisession, workers = 2L)

# read in the LAS Catalog -> create new LAS Catalog of classified points called GROUND 
GROUND = classify_ground(LCT, algorithm = csf(
  # in case there are steep slopes 
  sloop_smooth = TRUE,
  # minimum height threshold for canopy points
  class_threshold = 0.5,
  # larger resolution = coarser DTM
  cloth_resolution = 1,
  # rigidness of cloth - higher = more rigid
  rigidness = 3L,
  iterations = 500L,
  time_step = 1)) 

# alwasys stop the parallel processing!
future:::ClusterRegistry("stop")

#### CREATE A DEM FROM THE GROUND CLASSIFIED POINTS ####

# read in the LAS Catalog of GROUND classified points
GROUND = readLAScatalog(folder = paste0(dir, "\\output\\GROUND\\"))
opt_chunk_buffer(GROUND) = 15
# drop noisy points below the ground
opt_filter(GROUND) = "-drop_z_below -0.25"
opt_output_files(GROUND) = "" # no LAS output
opt_progress(GROUND) = TRUE

# create a 5m DEM from the GROUND points 
# using k-nn inverse distance weighting interpolation
DEM = grid_terrain(GROUND, res = 5, knnidw(), full_raster = TRUE)

#### NORMALIZATION ####
# normalize the point cloud to the ground classified points

# read the classified point cloud in with high-res parameters
GROUND = readLAScatalog(folder = paste0(dir, "\\output\\GROUND\\"))
# small buffer; we're just normalizing
opt_chunk_buffer(GROUND) = 1
opt_laz_compression(LCT) = TRUE
# output normalized point cloud to NORM folder
opt_output_files(GROUND) = paste0(dir, "\\output\\NORM\\{*}_NORM")
opt_progress(GROUND) = TRUE

# normalize the point cloud to the GROUND points
# or alternatively to the rasterized DEM - a bit faster but less "exact"
# NORM = normalize_height(LCT, DEM, knnidw(), na.rm = TRUE)
plan(multisession, workers = 2L)
NORM = normalize_height(GROUND, knnidw(), na.rm = TRUE)
future:::ClusterRegistry("stop")

## Create rasterized CHM from normalized point cloud
NORM = readLAScatalog(paste0(dir, "\\output\\NORM\\"))
# no buffer
opt_chunk_buffer(NORM) = 0
opt_laz_compression(NORM) = TRUE
# drop noisy points below the ground or floating in outer space
opt_filter(NORM) = "-drop_z_below -.25 -drop_z_above 30"
opt_output_files(NORM) = ""
opt_progress(NORM) = TRUE

# create a 1m CHM from the normalized point cloud using maximum canopy height
plan(multisession, workers = 2L)
CHM_max = pixel_metrics(NORM, res = 1, func = ~max(Z))
future:::ClusterRegistry("stop")