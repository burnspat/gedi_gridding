# Title: gedi_l2-l4a_s00_find_filter.R
# Author: Patrick Burns [pb463@nau.edu]
# Purpose: find GEDI sub-orbits and associated L2 and L4 granules which overlap a region of interest and filter by date
###############################################################################



#####
# Libraries
#####
source('/home/pb463/scripts/repos/gedi_gridding/R/functions.R')



#####
# Inputs
#####

# Get inputs passed from .sh file
args <- commandArgs(TRUE)

# GEDI L2 products Options are:
# L2A, L2B, L4A, L2A-L2B, L2A-L2B-L4A
prod=args[1]

# Name of study region
reg_name <- args[2]

# Geometry to use for clipping
clip <- args[3]

# Start date
sdate <- as.character(args[4])

# End date
edate <- as.character(args[5])

# Minimum DOY
min_doy = as.integer(args[6])

# Maximum DOY
max_doy = as.integer(args[7])

# Path to save files
save_file <- args[8]

# Paths for the L2A+B .h5 files
l2a_dir <- args[9] #'/projects/geode_data/share/lidar/GEDI02_A.002/0_orig/alldata/'
l2b_dir <- args[10] #'/projects/geode_data/share/lidar/GEDI02_B.002/0_orig/alldata/'
l4a_dir <- args[11] #'/projects/geode_data/share/lidar/GEDI04_A.002_1/0_orig/data/'



#####
# Processing
#####

# Handle roi input
# Need to convert to bounding box coordinates in LL Longitude, LL Latitude, UR Longitude, UR Latitude format
if(tools::file_ext(clip) == "shp"){
  clip_sp <- rgdal::readOGR(dsn = clip)
  clip_sp_tr <- sp::spTransform(x = clip_sp, CRSobj = sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  bbox <- paste0(clip_sp_tr@bbox[1],",", clip_sp_tr@bbox[2],",", clip_sp_tr@bbox[3],",", clip_sp_tr@bbox[4])
} else {
  inCoords <- as.numeric(unlist(strsplit(clip, ",")))
  bbox <- paste0(inCoords[2],",", inCoords[3],",", inCoords[4],",", inCoords[1])
}

# Run GEDI finder function to spatially query sub-orbit granules
if (grepl('L2A', prod, fixed = TRUE)){
  DT_2a <- data.table(gran_2a = paste0(l2a_dir, basename(unlist(gedi_finder(product = 'GEDI02_A.002', bbox = bbox)))))
  DT_2a <- DT_2a[, ext := tools::file_ext(gran_2a)]
  DT_2a <- DT_2a[ext == "h5"]
  DT_2a <- DT_2a[, key := substr(basename(DT_2a$gran_2a),10,32)]
}

if (grepl('L2B', prod, fixed = TRUE)){
  DT_2b <- data.table(gran_2b = paste0(l2b_dir, basename(unlist(gedi_finder(product = 'GEDI02_B.002', bbox = bbox)))))
  DT_2b <- DT_2b[, ext := tools::file_ext(gran_2b)]
  DT_2b <- DT_2b[ext == "h5"]
  DT_2b <- DT_2b[, key := substr(basename(DT_2b$gran_2b),10,32)]
}

if (grepl('L4A', prod, fixed = TRUE)){
  DT_4a <- data.table(gran_4a = paste0(l4a_dir, basename(unlist(gedi_finder(product = 'GEDI04_A.002', bbox = bbox)))))
  DT_4a <- DT_4a[, ext := tools::file_ext(gran_4a)]
  DT_4a <- DT_4a[ext == "h5"]
  DT_4a <- DT_4a[, key := substr(basename(DT_4a$gran_4a),10,32)]
}

# Join
if (prod == "L2A-L2B"){
  setkey(DT_2a, "key")
  setkey(DT_2b, "key")
  DT_s <- DT_2a[DT_2b, nomatch = 0]

} else if (prod == "L2A-L2B-L4A"){
  setkey(DT_2a, "key")
  setkey(DT_2b, "key")
  setkey(DT_4a, "key")
  DT_2ab <- DT_2a[DT_2b, nomatch = 0]
  DT_s <- DT_2ab[DT_4a, nomatch = 0]

} else if (prod == "L2A"){
  DT_s <- DT_2a

} else if (prod == "L2B"){
  DT_s <- DT_2b

} else if (prod == "L4A"){
  DT_s <- DT_4a
}

# Filter by date
sdate_d <- decimal_date(ymd(sdate))
edate_d <- decimal_date(ymd(edate))

DT_s <- DT_s[, ':='(origin = paste0(substr(key,1,4),"-01-01"), DOY = as.numeric(substr(key,5,7))) ]
DT_s <- DT_s[, date_d := decimal_date(as.Date(DOY-1, origin))]

DT_sf <- DT_s[date_d >= sdate_d & date_d <= edate_d & DOY >= min_doy & DOY <= max_doy]

cat("\n")
cat(nrow(DT_sf),"granules to process \n")

# Save filtered granules
if (prod == "L2A-L2B"){
  DT_sf <- DT_sf[, .(gran_2a, gran_2b)]

} else if (prod == "L2A-L2B-L4A"){
  DT_sf <- DT_sf[, .(gran_2a, gran_2b, gran_4a)]

} else if (prod == "L2A"){
  DT_sf <- DT_sf[, .(gran_2a)]

} else if (prod == "L2B"){
  DT_sf <- DT_sf[, .(gran_2b)]

} else if (prod == "L4A"){
  DT_sf <- DT_sf[, .(gran_4a)]
}

fwrite(x = DT_sf, file = save_file, row.names = FALSE, col.names = FALSE)
cat("Saved file list to", save_file, "\n")
cat("\n")
