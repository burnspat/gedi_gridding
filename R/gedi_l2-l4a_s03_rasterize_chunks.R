# Title: gedi_l2-l4a_s03_rasterize_chunks.R
# Author: Patrick Burns [pb463@nau.edu]
# Purpose: rasterize GEDI shots using several aggregation functions
###############################################################################



#####
# Libraries
source('/home/pb463/scripts/repos/gedi_gridding/R/functions.R')
cat("\n")
#####



#####
# 1. Inputs
#####

# Get inputs passed from .sh file
args <- commandArgs(TRUE)

# chunk id (e.g. G069W02N)
chunk_id <- args[1] 

# processed and combined GEDI data (including shots outside the chunk as buffer)
gedi_tab <- args[2]

# GEDI metrics
metrics <-unlist(strsplit(x = args[3], split = ","))

# date ranges
date_pairs <- unlist(strsplit(x = args[4], split = ","))

# base raster directory
base_raster_dir <- args[5]

# buffer distance
buff <- as.numeric(args[6])

# directory to save results
save_dir <- args[7]

# terra tmp directory (optional)
# terraOptions(tempdir = "/scratch/pb463/R_WD/terra_temp/")

# output file prefix
prefix <- "gediv002_"



#####
# 2. Processing
#####

# specify the output column order (n=158)
out_col_order <- c("shot_num", "orbit", "sub_orbit", "ease72_og", "beam", "beam_azim", "beam_elev", "delta_time", "date_dec", "date", "time", "doy", "sol_elev", "sol_azim", "lon_lm_a0", "lat_lm_a0", 
"elev_lm_a0", "elev_hr_a0", "elev_dem_srtm", "elev_dem_tdx", "elev_diff_dem", "stale_flag", "deg_flag", "surf_flag", "loc_out_umd", "l2a_selalg_a0", "l2a_qflag_a0", "l2a_qflag_a1", 
"l2a_qflag_a2", "l2a_qflag_a3", "l2a_qflag_a4", "l2a_qflag_a5", "l2a_qflag_a6", "l2_hqflag", "sens_a0", "sens_a1", "sens_a2", "sens_a3", "sens_a4", "sens_a5", "sens_a6", "sens_a10", 
"energy_total", "ls_treecov", "ls_waterp", "leafoff_flag", "leafoff_doy", "leafon_doy", "pft", "region", "urb_prop", "num_modes_a0", "rh_0_a0", "rh_5_a0", "rh_10_a0", "rh_15_a0", "rh_20_a0", 
"rh_25_a0", "rh_30_a0", "rh_35_a0", "rh_40_a0", "rh_45_a0", "rh_50_a0", "rh_55_a0", "rh_60_a0", "rh_65_a0", "rh_70_a0", "rh_75_a0", "rh_80_a0", "rh_85_a0", "rh_90_a0", "rh_92_a0", "rh_95_a0", 
"rh_98_a0", "rh_99_a0", "rh_100_a0", "rhvdr_b", "rhvdr_m", "rhvdr_t", "l2b_algrun_flag", "l2b_qflag_a0", "pgap_theta_a0", "pgap_theta_err_a0", "fhd_pai_1m_a0", "cover_a0", "cover_l1", "cover_l2", 
"cover_l3", "cover_l4", "cover_l5", "cover_l6", "cover_l7", "cover_l8", "cover_l9", "cover_l10", "cover_l11", "cover_l12", "cover_l13", "cover_l14", "cover_l15", "cover_l16", "cover_l17", "cover_l18", 
"cover_l19", "cover_l20", "cover_l21", "pai_a0", "pai_l1", "pai_l2", "pai_l3", "pai_l4", "pai_l5", "pai_l6", "pai_l7", "pai_l8", "pai_l9", "pai_l10", "pai_l11", "pai_l12", "pai_l13", "pai_l14", "pai_l15", 
"pai_l16", "pai_l17", "pai_l18", "pai_l19", "pai_l20", "pai_l21", "pavd_0_5", "pavd_5_10", "pavd_10_15", "pavd_15_20", "pavd_20_25", "pavd_25_30", "pavd_30_35", "pavd_35_40", "pavd_40_45", "pavd_45_50", 
"pavd_50_55", "pavd_55_60", "pavd_60_65", "pavd_65_70", "pavd_70_75", "pavd_75_80", "pavd_80_85", "pavd_85_90", "pavd_90_95", "pavd_95_100", "pavd_0_5_frac", "pavd_bot_frac", "pavd_top_frac", "pavd_max_h", 
"even_pai_1m_a0", "fhd_pavd_5m_a0", "even_pavd_5m_a0", "l4a_algrun_flag", "l4a_qflag_a0", "l4a_hqflag", "l4a_pred_strat", "agbd_a0", "agbd_se_a0", "agbd_pi_low_a0", "agbd_pi_upper_a0")

# get extents of the chunk
xs = substr(chunk_id, 5, 5)
if (xs == "W"){
  x_min <- (-1*as.numeric(substr(chunk_id, 2, 4)))
} else if(xs == "E"){
  x_min <- as.numeric(substr(chunk_id, 2, 4))
}
x_max <- as.numeric(x_min + 1)
x_min_buff <- as.numeric(x_min - buff)
x_max_buff <- as.numeric(x_max + buff)

ys = substr(chunk_id, 8, 8)
if (ys == "S"){
  y_min <- (-1*as.numeric(substr(chunk_id, 6, 7)))
} else if (ys == "N"){
  y_min <- as.numeric(substr(chunk_id, 6, 7))
}
y_max <- as.numeric(y_min + 1)
y_min_buff <- as.numeric(y_min - buff)
y_max_buff <- as.numeric(y_max + buff)

cat("----- Loading and filtering the combined shot table ----- \n")
# Load processed and combined GEDI shots and sort by date
DT_raw <- fread(gedi_tab)
setorder(DT_raw, cols = "date_dec")

# reorder columns
DT_raw <- setcolorder(DT_raw, out_col_order)

# TODO: remove this block ~~~~~
# need to fix two metrics
DT_raw <- DT_raw[, fhd_pavd_5m_a0 := NULL]
DT_raw <- DT_raw[, even_pavd_5m_a0 := NULL]

  # shannon's H and evenness of the 5 m PAVD profile
  shandiv_pavd5_func <- function(pavd_0_5, pavd_5_10, pavd_10_15, pavd_15_20,
                                 pavd_20_25, pavd_25_30, pavd_30_35, pavd_35_40,
                                 pavd_40_45, pavd_45_50, pavd_50_55, pavd_55_60,
                                 pavd_60_65, pavd_65_70, pavd_70_75, pavd_75_80,
                                 pavd_80_85, pavd_85_90, pavd_90_95, pavd_95_100){
    x <- c(pavd_0_5, pavd_5_10, pavd_10_15, pavd_15_20,
           pavd_20_25, pavd_25_30, pavd_30_35, pavd_35_40,
           pavd_40_45, pavd_45_50, pavd_50_55, pavd_55_60,
           pavd_60_65, pavd_65_70, pavd_70_75, pavd_75_80,
           pavd_80_85, pavd_85_90, pavd_90_95, pavd_95_100)

    # compute relative proportions
    x_pos <- x[!is.na(x) & x>0]
    x_sum <- sum(x_pos)
    x_rp <- x_pos/x_sum
    if (length(x_rp) > 1){ # need to have at least two bins for diversity
      h <- -1*(sum(x_rp*log(x_rp))) 
      return(h)
    } else if (length(x_rp) == 1){ # only 1 bin is 0 diversity
      return(0)
    } else {
      return(-9999)
    }
  }

  evenness_pavd5_func <- function(fhd_pavd_5m_a0, 
                                  pavd_0_5, pavd_5_10, pavd_10_15, pavd_15_20,
                                  pavd_20_25, pavd_25_30, pavd_30_35, pavd_35_40,
                                  pavd_40_45, pavd_45_50, pavd_50_55, pavd_55_60,
                                  pavd_60_65, pavd_65_70, pavd_70_75, pavd_75_80,
                                  pavd_80_85, pavd_85_90, pavd_90_95, pavd_95_100){
    
    x <- c(pavd_0_5, pavd_5_10, pavd_10_15, pavd_15_20,
           pavd_20_25, pavd_25_30, pavd_30_35, pavd_35_40,
           pavd_40_45, pavd_45_50, pavd_50_55, pavd_55_60,
           pavd_60_65, pavd_65_70, pavd_70_75, pavd_75_80,
           pavd_80_85, pavd_85_90, pavd_90_95, pavd_95_100)
    if (max(x, na.rm = TRUE) > 0){
      n_valid_bins <- as.numeric(max(which(x>0), na.rm = TRUE))
    } else {
      n_valid_bins <- 0
    }
    
    if (n_valid_bins > 1 && fhd_pavd_5m_a0 > -9999){
      return(fhd_pavd_5m_a0/log(n_valid_bins))
    } else if (n_valid_bins == 1 && fhd_pavd_5m_a0 > -9999){
      return(0)
    } else {
      return(-9999)
    }
  }

  DT_raw[, fhd_pavd_5m_a0:= fifelse(rh_100_a0 > 5, 
                                    mapply(shandiv_pavd5_func, 
                                    pavd_0_5, pavd_5_10, pavd_10_15, pavd_15_20,
                                    pavd_20_25, pavd_25_30, pavd_30_35, pavd_35_40,
                                    pavd_40_45, pavd_45_50, pavd_50_55, pavd_55_60,
                                    pavd_60_65, pavd_65_70, pavd_70_75, pavd_75_80,
                                    pavd_80_85, pavd_85_90, pavd_90_95, pavd_95_100),
                                    -9999)]
  
  DT_raw[, even_pavd_5m_a0 := fifelse(rh_100_a0 > 5,
                                      mapply(evenness_pavd5_func, 
                                      fhd_pavd_5m_a0,
                                      pavd_0_5, pavd_5_10, pavd_10_15, pavd_15_20,
                                      pavd_20_25, pavd_25_30, pavd_30_35, pavd_35_40,
                                      pavd_40_45, pavd_45_50, pavd_50_55, pavd_55_60,
                                      pavd_60_65, pavd_65_70, pavd_70_75, pavd_75_80,
                                      pavd_80_85, pavd_85_90, pavd_90_95, pavd_95_100), 
                                      -9999)]
# ~~~~~

# subset the DT by buffered bounding box
DT_s <- DT_raw[lon_lm_a0 >= x_min_buff & lon_lm_a0 <= x_max_buff 
               & lat_lm_a0 >= y_min_buff & lat_lm_a0 <= y_max_buff]
cat("Subset GEDI DT to buffered bounding box:", x_min_buff, x_max_buff, y_min_buff, y_max_buff, "\n")
nrow_orig <- nrow(DT_s)
cat(nrow_orig, "shots \n")
rm(DT_raw)

# omit any rows which have NA values for coordinates, date, or elev_lm_a0
DT_ll <- na.omit(DT_s, cols=c("lon_lm_a0", "lat_lm_a0", "lon_lm_a0_6933", "lat_lm_a0_6933", "elev_lm_a0",
                              "date_dec", "orbit"))
nrow_ll_omit <- as.numeric(nrow(DT_ll) - nrow_orig)
cat("Omitted", nrow_ll_omit, "shots due to missing coordinates and/or decimal date \n")
rm(DT_s)

# remove shots/granules identified as local outliers
# DT_g will be used for subsequent aggregation of ground metrics
DT_g <- DT_ll[loc_out_umd == 0] 
nrow_bg <- as.numeric(nrow(DT_ll) - nrow(DT_g))
cat("Omitted", nrow_bg, "shots associated with outlier granules \n")

# for veg. metrics use the vegetation high quality flag
# DT_v will be used for subsequent aggregation of veg. metrics
DT_v <- DT_ll[l2_hqflag == 1]
nrow_lq <- as.numeric(nrow(DT_ll)-nrow(DT_v))
cat(nrow_lq, "shots are NOT high quality \n")
rm(DT_ll)

# function to compute summary statistics of numeric GEDI metrics
gedi_dt_summ <- function(dt, grid_id){
  dt_classes <- sapply(dt, class)
  numeric_cols <- unique(names(dt_classes[dt_classes == "numeric" | dt_classes == "integer"]))
  dt_summ <- data.table(grid_id = rep(x = grid_id, length(numeric_cols)),
                        field_name = numeric_cols,
                        n_all = sapply(dt[,..numeric_cols], function(x) length(x)),
                        n_notna = sapply(dt[,..numeric_cols], function(x) sum(!is.na(x))),
                        n_notnafin = sapply(dt[,..numeric_cols], function(x) sum(!is.finite(x) & !is.na(x))),
                        n_valid = sapply(dt[,..numeric_cols], function(x) sum(!is.na(x) & x!=-9999 & is.finite(x))),
                        n_n9999 = sapply(dt[,..numeric_cols], function(x) sum(x==-9999 & is.finite(x), na.rm=TRUE)),
                        n_255 = sapply(dt[,..numeric_cols], function(x) sum(x==255 & is.finite(x), na.rm=TRUE)),
                        n_32767 = sapply(dt[,..numeric_cols], function(x) sum(x==32767 & is.finite(x), na.rm=TRUE)),
                        n_neg = sapply(dt[,..numeric_cols], function(x) sum(x<0 & is.finite(x), na.rm = TRUE)),
                        n_neg_valid = sapply(dt[,..numeric_cols], function(x) sum(x<0 & x!=-9999 & is.finite(x), na.rm = TRUE)),
                        min_nna = sapply(dt[,..numeric_cols], function(x) min(x[is.finite(x)], na.rm=TRUE)),
                        median_nna = sapply(dt[,..numeric_cols], function(x) median(x[is.finite(x)], na.rm=TRUE)),
                        max_nna = sapply(dt[,..numeric_cols], function(x) max(x[is.finite(x)], na.rm=TRUE)),
                        min_valid = sapply(dt[,..numeric_cols], function(x) min(x[x!=-9999 & is.finite(x)], na.rm=TRUE)),
                        mean_valid = sapply(dt[,..numeric_cols], function(x) mean(x[x!=-9999 & is.finite(x)], na.rm=TRUE)),
                        median_valid = sapply(dt[,..numeric_cols], function(x) median(x[x!=-9999 & is.finite(x)], na.rm=TRUE)),
                        max_valid = sapply(dt[,..numeric_cols], function(x) max(x[x!=-9999 & is.finite(x)], na.rm=TRUE)),
                        sd_valid = sapply(dt[,..numeric_cols], function(x) sd(x[x!=-9999 & is.finite(x)], na.rm=TRUE)),
                        p05_valid = sapply(dt[,..numeric_cols], function(x) quantile(x = x[x!=-9999 & is.finite(x)], probs = c(0.05), names = FALSE, na.rm=TRUE)),
                        p25_valid = sapply(dt[,..numeric_cols], function(x) quantile(x = x[x!=-9999 & is.finite(x)], probs = c(0.25), names = FALSE, na.rm=TRUE)),
                        p75_valid = sapply(dt[,..numeric_cols], function(x) quantile(x = x[x!=-9999 & is.finite(x)], probs = c(0.75), names = FALSE, na.rm=TRUE)),
                        p95_valid = sapply(dt[,..numeric_cols], function(x) quantile(x = x[x!=-9999 & is.finite(x)], probs = c(0.95), names = FALSE, na.rm=TRUE))
  )
  dt_summ <- dt_summ[ ,iqr_valid := (p75_valid - p25_valid)]
  return(dt_summ)
}

# compute and save summary stats for the corresponding all ground csv
DT_gs <- DT_g[lon_lm_a0 >= x_min & lon_lm_a0 < x_max & lat_lm_a0 >= y_min & lat_lm_a0 < y_max]
DT_gs <- DT_gs[,..out_col_order]
# Change a few of the column nodata values
DT_gs <- DT_gs[leafoff_flag == 255, leafoff_flag:= -9999]
DT_gs <- DT_gs[leafoff_doy == 32767, leafoff_doy:= -9999]
DT_gs <- DT_gs[leafon_doy == 32767, leafon_doy:= -9999]
DT_gs <- DT_gs[ls_treecov == 255, ls_treecov:= -9999]
DT_gs <- DT_gs[pft == 255, pft:= -9999]
DT_gs <- DT_gs[ls_waterp == 255, ls_waterp:= -9999]
DT_gs <- DT_gs[urb_prop == 255, urb_prop:= -9999]

DT_g_summ <- gedi_dt_summ(DT_gs, chunk_id)
fwrite(x = DT_g_summ, file = paste0(save_dir, prefix, "l2l4a_ga_summ_", tolower(chunk_id), ".csv"))
cat("Computed summary stats (all lower quality ground shots) for this chunk.\n")
rm(DT_g_summ)
rm(DT_gs)

# save high quality veg. shots
if (nrow(DT_v) > 0){
  save_file <- paste0(save_dir, prefix, "l2l4a_va_", tolower(chunk_id), ".csv")
  if (!file.exists(save_file)){
    # Save the table only to the extent of the chunk
    DT_c <- DT_v[lon_lm_a0 >= x_min & lon_lm_a0 < x_max & lat_lm_a0 >= y_min & lat_lm_a0 < y_max]
    DT_c <- DT_c[,..out_col_order]

    # Change a few of the column nodata values
    DT_c <- DT_c[leafoff_flag == 255, leafoff_flag:= -9999]
    DT_c <- DT_c[leafoff_doy == 32767, leafoff_doy:= -9999]
    DT_c <- DT_c[leafon_doy == 32767, leafon_doy:= -9999]
    DT_c <- DT_c[ls_treecov == 255, ls_treecov:= -9999]
    DT_c <- DT_c[pft == 255, pft:= -9999]
    DT_c <- DT_c[ls_waterp == 255, ls_waterp:= -9999]
    DT_c <- DT_c[urb_prop == 255, urb_prop:= -9999]

    fwrite(x = DT_c, file = save_file)
    cat("Saved", nrow(DT_c), "shots corresponding only to this chunk (no buffer) to", save_file, "\n")
    
    # compute and save summary stats for the corresponding csv
    DT_c_summ <- gedi_dt_summ(DT_c, chunk_id)
    fwrite(x = DT_c_summ, file = paste0(save_dir, prefix, "l2l4a_va_summ_", tolower(chunk_id), ".csv"))
    cat("Computed summary stats (all high quality veg. shots) for this chunk.\n")

    rm(DT_c)
    rm(DT_c_summ)
    cat("\n\n")
  } else{
    cat("File of shots corresponding only to this chunk (no buffer) already exists \n")
  }
} else {
  cat("0 shots for veg. rasterization. Quitting... \n")
  cat("\n")
}


# load the base rasters
base_rasters <- list.files(path = base_raster_dir, pattern = "*.tif", full.names = TRUE, recursive = FALSE)
base_epsg <- as.character(crs(rast(base_rasters[1]), describe=TRUE)[3])
cat("The base raster EPSG code is", base_epsg, "\n")
cat("\n")

if (base_epsg == "4326"){
  # identify GEDI coords
  coords <- c("lon_lm_a0", "lat_lm_a0")

  # make a 30 m raster for removing redundant shots
  fine_res <- 0.00026949335
  r_sac <- terra::rast(xmin = (x_min_buff - fine_res), xmax = (x_max_buff + fine_res), ymin = (y_min_buff - fine_res), ymax = (y_max_buff + fine_res),
                       crs = "epsg:4326", resolution = fine_res)

} else if (base_epsg == "6933"){
  coords <- c("lon_lm_a0_6933", "lat_lm_a0_6933")
  
  # Transform xmin, xmax, ymin, ymax
  sf_pt_ll <- st_geometry(st_point(c(x_min, y_min)))
  st_crs(sf_pt_ll) <- 4326
  coords_ll_6933 <- st_transform(x = sf_pt_ll, crs = 6933)
  sf_pt_ur <- st_geometry(st_point(c(x_min+1, y_min+1)))
  st_crs(sf_pt_ur) <- 4326
  coords_ur_6933 <- st_transform(x = sf_pt_ur, crs = 6933)
  x_min <- st_coordinates(coords_ll_6933)[1]
  y_min <- st_coordinates(coords_ll_6933)[2]
  x_max <- st_coordinates(coords_ur_6933)[1]
  y_max <- st_coordinates(coords_ur_6933)[2]

  # Transform buffered xmin, xmax, ymin, ymax
  sf_pt_ll_buff <- st_geometry(st_point(c(x_min_buff, y_min_buff)))
  st_crs(sf_pt_ll_buff) <- 4326
  coords_ll_6933_buff <- st_transform(x = sf_pt_ll_buff, crs = 6933)
  sf_pt_ur_buff <- st_geometry(st_point(c(x_max_buff, y_max_buff)))
  st_crs(sf_pt_ur_buff) <- 4326
  coords_ur_6933_buff <- st_transform(x = sf_pt_ur_buff, crs = 6933)
  x_min_buff <- st_coordinates(coords_ll_6933_buff)[1]
  y_min_buff <- st_coordinates(coords_ll_6933_buff)[2]
  x_max_buff <- st_coordinates(coords_ur_6933_buff)[1]
  y_max_buff <- st_coordinates(coords_ur_6933_buff)[2]

  # make a 30 m raster for removing redundant shots
  fine_res <- 30
  Round <- function(x,y) {
    if((y - x %% y) <= x %% y) { x + (y - x %% y)} 
    else { x - (x %% y)}
  }
  x_min_buff_r <- as.numeric(Round((x_min_buff - fine_res), 30))
  y_min_buff_r <- as.numeric(Round((y_min_buff - fine_res), 30))
  x_max_buff_r <- as.numeric(Round((x_max_buff + fine_res), 30))
  y_max_buff_r <- as.numeric(Round((y_max_buff + fine_res), 30))
  r_sac <- terra::rast(xmin = x_min_buff_r , xmax = x_max_buff_r, ymin = y_min_buff_r, ymax = y_max_buff_r,
                       crs = "epsg:6933", resolution = fine_res)
  # values(r_sac) <- runif(ncell(r_sac))
  # writeRaster(x = r_sac, filename = paste0(save_dir, prefix, "rsac.tif"), 
  #             gdal=c("COMPRESS=LZW", "TILED=YES"), NAflag = -9999, overwrite = TRUE)
} else {
  cat("GEDI shots not available in this CRS (", base_epsg, ")\n")
}

# format dates
dates_list <- list()
for (d in 1: length(date_pairs)){
  if ((d %% 2) == 1){
    date_time <- paste0(date_pairs[d], " 00:00:00")
    date_dec <- decimal_date(ymd_hms(date_time))
  } else if ((d %% 2) == 0){
    date_time <- paste0(date_pairs[d], " 23:59:59")
    date_dec <- decimal_date(ymd_hms(date_time))
  }
  dates_list[[d]] <- date_dec
}


# custom functions for raster aggregation
# interquartile range
iqr_func <- function(x){
  x <- x[!is.na(x) & is.finite(x)]
  percs <- quantile(x, probs = c(0.25, 0.75), names = FALSE, na.rm = TRUE)
  iqr <- as.numeric(percs[2]-percs[1])
  return(iqr)
}

# 95th percentile
p95_func <- function(x){
  x <- x[!is.na(x) & is.finite(x)]
  p95 <- quantile(x, probs = c(0.95), names = FALSE, na.rm = TRUE)
  return(p95)
}

# shannons H
# bins for metrics of interest
# using conservative maximum values for metrics which are unbounded
bins_min_dict <- c("date_dec" = 2019.0, "elev_lm_a0" = -200, "num_modes_a0" = 1, "sens_a0" = 0.9, "rh_50_a0" = -99, "rh_95_a0" = -99, "rh_98_a0" = -99,
                   "pai_a0" = 0, "fhd_pai_1m_a0" = 0, "cover_a0" = 0, "agbd_a0" = 0, "pavd_0_5_frac" = 0, "pavd_max_h" = 0,
                   "rhvdr_b" = 0, "rhvdr_m" = 0, "rhvdr_t" = 0, "pavd_bot_frac" = 0, "pavd_top_frac" = 0,
                   "fhd_pavd_5m_a0" = 0, "even_pai_1m_a0" = 0, "even_pavd_5m_a0" = 0, 
                   "pavd_0_5" = 0, "pavd_5_10" = 0, "pavd_10_15" = 0, "pavd_15_20" = 0, "pavd_20_25" = 0, "pavd_25_30" = 0,
                   "pavd_30_35" = 0, "pavd_35_40" = 0, "pavd_40_45" = 0, "pavd_45_50" = 0, "pavd_50_55" = 0, "pavd_55_60" = 0,
                   "pavd_60_65" = 0, "pavd_65_70" = 0, "pavd_70_75" = 0, "pavd_75_80" = 0)
bins_max_dict <- c("date_dec" = 2023.5, "elev_lm_a0" = 9000, "num_modes_a0" = 20, "sens_a0" = 1.0, "rh_50_a0" = 120, "rh_95_a0" = 120, "rh_98_a0" = 120,
                   "pai_a0" = 12, "fhd_pai_1m_a0" = 4, "cover_a0" = 1.0, "agbd_a0" = 8000, "pavd_0_5_frac" = 1, "pavd_max_h" = 100,
                   "rhvdr_b" = 1, "rhvdr_m" = 1, "rhvdr_t" = 1, "pavd_bot_frac" = 1, "pavd_top_frac" = 1,
                   "fhd_pavd_5m_a0" = 4, "even_pai_1m_a0" = 2, "even_pavd_5m_a0" = 1, 
                   "pavd_0_5" = 1.0, "pavd_5_10" = 1.0, "pavd_10_15" = 1.0, "pavd_15_20" = 1.0, "pavd_20_25" = 1.0, "pavd_25_30" = 1.0,
                   "pavd_30_35" = 1.0, "pavd_35_40" = 1.0, "pavd_40_45" = 1.0, "pavd_45_50" = 1.0, "pavd_50_55" = 1.0, "pavd_55_60" = 1.0,
                   "pavd_60_65" = 1.0, "pavd_65_70" = 1.0, "pavd_70_75" = 1.0, "pavd_75_80" = 1.0)
bins_step_dict <- c("date_dec" = 0.08333333, "elev_lm_a0" = 200, "num_modes_a0" = 1, "sens_a0" = 0.005, "rh_50_a0" = 1.5, "rh_95_a0" = 3, "rh_98_a0" = 3,
                    "pai_a0" = 0.25, "fhd_pai_1m_a0" = 0.1, "cover_a0" = 0.05, "agbd_a0" = 20, "pavd_0_5_frac" = 0.025, "pavd_max_h" = 5,
                    "rhvdr_b" = 0.025, "rhvdr_m" = 0.025, "rhvdr_t" = 0.025, "pavd_bot_frac" = 0.025, "pavd_top_frac" = 0.025,
                    "fhd_pavd_5m_a0" = 0.1, "even_pai_1m_a0" = 0.01, "even_pavd_5m_a0" = 0.05, 
                    "pavd_0_5" = 0.01, "pavd_5_10" = 0.01, "pavd_10_15" = 0.01, "pavd_15_20" = 0.01, "pavd_20_25" = 0.01, "pavd_25_30" = 0.01,
                    "pavd_30_35" = 0.01, "pavd_35_40" = 0.01, "pavd_40_45" = 0.01, "pavd_45_50" = 0.01, "pavd_50_55" = 0.01, "pavd_55_60" = 0.01,
                    "pavd_60_65" = 0.01, "pavd_65_70" = 0.01, "pavd_70_75" = 0.01, "pavd_75_80" = 0.01)

shandiv_func <- function(x, bins){
  x <- x[!is.na(x) & is.finite(x)] #& x>=0
  if (length(x)>1){
    hist <- hist(x = x, breaks = bins, plot = FALSE) 
    p_x <- hist$counts[hist$counts>0]/sum(hist$counts)

    if (length(p_x) > 1){ # need to have at least two bins for diversity
      h <- -1*(sum(p_x*log(p_x))) 
      return(h)
    } else if (length(p_x) == 1){ # only 1 bin is 0 diversity
      return(0)
    } else {
      return(NA)
    }
  } else {
    return(NA)
  }
}

# number of shots excluding NA
count_func <- function(x){
  shot_count <- length(x = na.omit(x))
  return(shot_count)
}

# number of unique values excluding NA
n_uniq_func <- function(x){
  n_uniq <- uniqueN(x = x, na.rm = TRUE)
  return(n_uniq)
}

# Average Nearest Neighbor Index (NNI)
# adapted from R spatialEco (https://github.com/jeffreyevans/spatialEco/blob/master/R/nni.R)
nni_func <- function(x, y, area){
  ppp <- ppp(x = x, y = y, unitname=c("metre","metres"), check = FALSE)
  obsMeanDist <- sum(spatstat.geom::nndist(ppp))/ppp$n
  expMeanDist <- 0.5 * sqrt(area / ppp$n)
  nni <- as.numeric(obsMeanDist/expMeanDist)
  if (is.finite(nni)){
    return(nni)
  } else {
    return(-9999)
  }
}

# bootstrap standard error of the mean
# first try - doesn't guarantee unique bootstrap samples, esp. for cells with few shots
# mean_boot_se <- function(x){
#   x <- x[!is.na(x) & is.finite(x)]
#   if (length(x)>=10){
#     ss <- round(length(x)*0.7)
#     boots <- seq(1,100)
#     for (b in boots) {
#       set.seed(b)
#       bs <- sample(x = x, size = ss, replace=FALSE)
#       boots[b] <- mean(bs)
#     }
#     boot_se <- sd(boots)
#     return(boot_se)
#   } else {
#     return(NA)
#   }
# }

mean_boot_se <- function(x){
  x <- x[!is.na(x) & is.finite(x)]
  if (length(x)>=10){
    # specify the proportion of shots to subsample
    ss <- round(length(x)*0.7)

    # specify a reasonable number of bootstraps considering the number of shots
    # need to use more bootstraps when there are fewer shots in order to ensure at least 100 unique boots
    if (length(x)>12){
      n_boots <- 150
    } else {
      n_boots <- 400
    }

    # make list of random sample indices 
    bsi_list <- list()
    for (b in 1:n_boots){
      set.seed(b)
      bsi <- sort(sample(x = seq(1, length(x)), size = ss, replace = FALSE))
      bsi_list[[b]] <- bsi
    }

    # take the first 100 unique index lists
    if (length(unique(bsi_list)) >= 100){
      bsi_list_u <- unique(bsi_list)[1:100]
    } else {
      # just in case there aren't 100 unique index lists, use all unique that are available
      bsi_list_u <- unique(bsi_list)
    }

    #  get the values associated with each index list and make a list of vectors and calculate the SE
    se <- sd(unlist(lapply(bsi_list_u, function(i){return(mean(x[i]))})))
    return(se) 
  } else {
    return(NA)
  }
}


# add a variation of the agbd_a0 metric which uses the L4a quality flag
if ("agbd_a0" %in% metrics){
  metrics <- c(metrics, "agbd_a0_qf")
}

# statistics dictionary to use for rasterization
stats_dict <- c("mean"="mean", "meanbse"=mean_boot_se, "med"="median", "sd"="sd", "iqr"=iqr_func, "p95"=p95_func, "shan"=shandiv_func, "countf"=count_func)

# rasterize by spatial resolution (b), dates(d), metrics (m), and stats (s)
for (b in 1:length(base_rasters)){
  base <- rast(base_rasters[b])

  # get the spatial resolution in meters
  if (base_epsg == "4326"){
    res <- round(xres(base)*111320)
    cell_area <- (xres(base)*111320)^2 #TODO: this is a rough approximation (but we shouldn't be calculating areas for geographic coordinates anyway)
  } else if (base_epsg == "6933"){
    res <- round(xres(base))
    cell_area <- (xres(base))^2
  }
  cat("----- Begin rasterization at", res, "m spatial resolution ----- \n")
  
  # get the extent for cropping (+2 pixels on all sides)
  x_buffer <- 1*as.numeric(xres(base))
  y_buffer <- 1*as.numeric(yres(base))
  x_min_wind <- (x_min - x_buffer)
  x_max_wind <- (x_max + x_buffer)
  y_min_wind <- (y_min - y_buffer)
  y_max_wind <- (y_max + y_buffer)
  ext_crop <- ext(x_min_wind, x_max_wind, y_min_wind, y_max_wind)

  # Subset the DT by start and stop date pairs
  di <- 1
  for (d in 1:(length(dates_list)/2)){
    date_dec_start <- as.numeric(dates_list[[di]])
    date_dec_start_lab <- gsub(pattern = "-", replacement = "", x = as.character(as_date(date_decimal(decimal = date_dec_start))))
    date_dec_end <- as.numeric(dates_list[[di+1]])
    date_dec_end_lab <- gsub(pattern = "-", replacement = "", x = as.character(as_date(date_decimal(decimal = date_dec_end))))
    
    cat(" ~~~ Begin Date Range: ", date_dec_start_lab, "to", date_dec_end_lab, "~~~\n")
    
    # temporally filter tables for veg and ground metrics 
    # veg metrics
    DT_tv <- DT_v[date_dec >= date_dec_start & date_dec <= date_dec_end]
    # ground
    DT_tg <- DT_g[date_dec >= date_dec_start & date_dec <= date_dec_end]

    # extract 30 m resolution raster cell values 
    DT_tv <- DT_tv[ ,fine_cell:=cellFromXY(object = r_sac, xy = .SD), .SDcols = coords]
    DT_tg <- DT_tg[ ,fine_cell:=cellFromXY(object = r_sac, xy = .SD), .SDcols = coords]

    # get the first shot per 30 m cell
    DT_tvf <- DT_tv[, .SD[1], by = fine_cell]
    DT_tgf <- DT_tg[, .SD[1], by = fine_cell]

    if (nrow(DT_tg) > 0){
      cat(">> Working on count metrics - there are", nrow(DT_tg), "quality ground shots for count rasterization \n")
      cat(">>>> Calculating GEDI ground shot count \n")
      agg_shotct_g <- rasterize(x = as.matrix(DT_tg[, ..coords]), values = DT_tg[, date_dec], y = base, fun = count_func)
      names(agg_shotct_g) <- paste0("shots_count")
      shot_ck_g <- ifelse(uniqueN(x = values(agg_shotct_g), na.rm = TRUE) > 0, TRUE, FALSE)
      shotct_mask_g <- agg_shotct_g
      shotct_mask_g[shotct_mask_g <= 1] <- NA
      
      cat(">>>> Calculating GEDI ground shot NNI \n")
      #TODO: verirfy this data.table calc and load spatstat in gedi_funcs.R
      DT_tg_nni <- DT_tg[, cell_n:=cellFromXY(object = base,  xy = .SD), .SDcols = coords]
      DT_tg_nni <- DT_tg_nni[, nni:= mapply(nni_func, .SD[,1], .SD[,2], cell_area), by = cell_n, .SDcols = coords]
      DT_tg_nni <- DT_tg_nni[nni > -9999]
      agg_shotnni_g <- rasterize(x = as.matrix(DT_tg_nni[, ..coords]), values = DT_tg_nni[, nni], y = base, fun = "max")
      names(agg_shotnni_g ) <- paste0("shots_nni")
      rm(DT_tg_nni)

      cat(">>>> Calculating GEDI ground orbit count \n")
      agg_orbitct_g <- rasterize(x = as.matrix(DT_tg[, ..coords]), values = DT_tg[, orbit], y = base, fun = n_uniq_func)
      names(agg_orbitct_g) <- paste0("orbits_uniq")
      orbit_ck_g <- ifelse(uniqueN(x = values(agg_orbitct_g), na.rm = TRUE) > 0, TRUE, FALSE)

      cat(">>>> Calculating GEDI ground track count \n")
      DT_tg <- DT_tg[, track := paste0(orbit,"-",beam)]
      agg_trackct_g <- rasterize(x = as.matrix(DT_tg[, ..coords]), values = DT_tg[, track], y = base, fun = n_uniq_func)
      names(agg_trackct_g) <- paste0("tracks_uniq")
      track_ck_g <- ifelse(uniqueN(x = values(agg_trackct_g), na.rm = TRUE) > 0, TRUE, FALSE)
      
      # check that there is some data before writing
      if (any(shot_ck_g, orbit_ck_g, track_ck_g)){
        # stack the count rasters
        agg_counts_g <- rast(list(agg_shotct_g, agg_orbitct_g, agg_trackct_g, agg_shotnni_g))
        cat(">> Saving all ground count raster \n")
        writeRaster(x = crop(agg_counts_g, ext_crop), filename = paste0(save_dir, prefix, "counts_ga_", date_dec_start_lab, "_", date_dec_end_lab, "_", as.character(res), "m.tif"), 
                    gdal=c("COMPRESS=LZW", "TILED=YES"), NAflag = -9999, overwrite = TRUE)
        rm(agg_orbitct_g, agg_trackct_g, agg_counts_g, agg_shotnni_g)
      } else {
        cat(">> NOT saving all ground count raster because there are 0 valid values \n")
      }
      cat("\n")

      cat(">> Working on first elev_lm_a0 stats - there are", nrow(DT_tgf), "shots for rasterization \n")
      m_rast_list <- list()
      m_ck_list <- list()    
      
      for (s in seq(1, length(stats_dict))){
        stat_name <- as.character(names(stats_dict[s]))
        cat(">>>>", stat_name, "\n")
        if (stat_name == "shan"){
          # define bins for Shannon's H
          bins_min <- as.numeric(bins_min_dict["elev_lm_a0"])
          bins_max <- as.numeric(bins_max_dict["elev_lm_a0"])
          bins_step <- as.numeric(bins_step_dict["elev_lm_a0"])
          bins <- seq(bins_min, bins_max, bins_step)
          
          agg_stat <- rasterize(x = as.matrix(DT_tgf[, ..coords]), values = DT_tgf[,elev_lm_a0], 
                                y = base, fun = stats_dict[[s]], bins = bins)
        } else {
          agg_stat <- rasterize(x = as.matrix(DT_tgf[, ..coords]), values = DT_tgf[,elev_lm_a0], 
                                y = base, fun = stats_dict[[s]])
        }

        names(agg_stat) <- paste0("elev-lm-a0_", as.character(names(stats_dict[s])))
        m_rast_list[[s]] <- agg_stat
          
        # check for valid values 
        m_ck <-ifelse(uniqueN(x = values(agg_stat), na.rm = TRUE) > 0, TRUE, FALSE)
        m_ck_list[[s]] <- m_ck
      } # ends the stats loop

      # check that at least 1 stat raster has valid values
      if (any(unlist(m_ck_list))){ 
        # stack all of the statistics for the metric of interest and export
        agg_stats_gf <- rast(m_rast_list)
        shotctgf_mask <- agg_stats_gf["elev-lm-a0_countf"]
        shotctgf_mask[shotctgf_mask <= 1] <- NA
        
        # mask out pixels that have 1 shot or less
        agg_stats_gf_masked = mask(agg_stats_gf, shotctgf_mask)
        
        cat(">> Saving first elev_lm_a0 stats raster \n")
        writeRaster(x = crop(agg_stats_gf_masked, ext_crop), filename = paste0(save_dir, prefix, "elev-lm-a0_gf_", date_dec_start_lab, "_", date_dec_end_lab, "_", as.character(res), "m.tif"), 
                    gdal=c("COMPRESS=LZW", "TILED=YES"), NAflag = -9999, overwrite = TRUE)
        rm(agg_stat, m_rast_list, agg_stats_gf, agg_stats_gf_masked)
      } else {
        cat(">>> NOT saving first elev_lm_a0 stats raster because there are 0 valid values \n") 
      }
    } else { # if there are no shots within this time period
      cat("  there are 0 elev_lm_a0 shots to rasterize during this time period \n")
    }
    rm(DT_tg)
    rm(DT_tgf)
    cat("\n")

    # check that there are veg. shots for this time window and extent
    if (nrow(DT_tv) > 0){
      cat(">> Working on vegetation count metrics - there are", nrow(DT_tv), "high quality veg. shots for rasterization \n")
      cat(">>>> Calculating GEDI veg. shot count \n")
      agg_shotct <- rasterize(x = as.matrix(DT_tv[, ..coords]), values = DT_tv[, date_dec], y = base, fun = count_func)
      names(agg_shotct) <- paste0("shots_count")
      shot_ck <- ifelse(uniqueN(x = values(agg_shotct), na.rm = TRUE) > 0, TRUE, FALSE)
      # make a shot count mask for masking metric stats
      shotct_mask <- agg_shotct
      shotct_mask[shotct_mask <= 1] <- NA

      cat(">>>> Calculating GEDI veg. shot NNI \n")
      #TODO: verirfy this data.table calc and load spatstat in gedi_funcs.R
      DT_tv_nni <- DT_tv[, cell_n:=cellFromXY(object = base,  xy = .SD), .SDcols = coords]
      DT_tv_nni <- DT_tv_nni[, nni:= mapply(nni_func, .SD[,1], .SD[,2], cell_area), by = cell_n, .SDcols = coords]
      DT_tv_nni <- DT_tv_nni[nni > -9999]
      agg_shotnni <- rasterize(x = as.matrix(DT_tv_nni[, ..coords]), values = DT_tv_nni[, nni], y = base, fun = "max")
      names(agg_shotnni ) <- paste0("shots_nni")
      rm(DT_tv_nni)

      cat(">>>> Calculating GEDI veg. orbit count \n")
      agg_orbitct <- rasterize(x = as.matrix(DT_tv[, ..coords]), values = DT_tv[, orbit], y = base, fun = n_uniq_func)
      names(agg_orbitct) <- paste0("orbits_uniq")
      orbit_ck <- ifelse(uniqueN(x = values(agg_orbitct), na.rm = TRUE) > 0, TRUE, FALSE)
      
      cat(">>>> Calculating GEDI veg. track count \n")
      DT_tv <- DT_tv[, track := paste0(orbit,"-",beam)]
      agg_trackct <- rasterize(x = as.matrix(DT_tv[, ..coords]), values = DT_tv[, track], y = base, fun = n_uniq_func)
      names(agg_trackct) <- paste0("tracks_uniq")
      track_ck <- ifelse(uniqueN(x = values(agg_trackct), na.rm = TRUE) > 0, TRUE, FALSE)

      if (any(shot_ck, orbit_ck, track_ck)){
        # stack the count rasters
        agg_counts <- rast(list(agg_shotct, agg_orbitct, agg_trackct, agg_shotnni))
        cat(">> Saving all vegetation count raster \n")
        writeRaster(x = crop(agg_counts, ext_crop), filename = paste0(save_dir, prefix, "counts_va_", date_dec_start_lab, "_", date_dec_end_lab, "_", as.character(res), "m.tif"), 
                    gdal=c("COMPRESS=LZW", "TILED=YES"), NAflag = -9999, overwrite = TRUE)
        rm(agg_orbitct, agg_trackct, agg_counts, agg_shotnni)
      } else {
        cat(">> NOT saving all vegetation count raster because there are 0 valid values \n")
      }
      cat("\n")

      # stats for each metric
      for (m in metrics){
        m_rast_list <- list()
        m_ck_list <- list()
        m_lab <- gsub(pattern = "_", replacement = "-", x = m)

        if (m == "agbd_a0_qf"){
          # apply the L4a quality flag
          DT_r <- DT_tvf[l4a_hqflag == 1]
          m <- "agbd_a0"
        } else {
          # use already applied l2a_hqflag shots for veg. metrics
          DT_r <- DT_tvf
        }

        # exclude no data values and infinite values
        DT_r <- DT_r[get(m) > -9999]
        DT_r <- DT_r[is.finite(get(m))]
        DT_r <- na.omit(DT_r, cols = c(m))

        cat(">> Calculating GEDI metric", m_lab, "stats \n")
        
        for (s in seq(1, length(stats_dict))){
          stat_name <- as.character(names(stats_dict[s]))
          cat(">>>>", stat_name, "\n")
          if (stat_name == "shan"){
            # define bins for Shannon's H
            bins_min <- as.numeric(bins_min_dict[m])
            bins_max <- as.numeric(bins_max_dict[m])
            bins_step <- as.numeric(bins_step_dict[m])
            bins <- seq(bins_min, bins_max, bins_step)

            # check that all values are in the range of the bins
            min_val <- min(DT_r[,get(m)])
            max_val <- max(DT_r[,get(m)])
            n_vals_out <- nrow(DT_r[get(m)<min(bins_min) | get(m)>max(bins_max)])
            if (n_vals_out == 0){
              agg_stat <- rasterize(x = as.matrix(DT_r[, ..coords]), values = DT_r[,get(m)], y = base, fun = stats_dict[[s]], bins = bins)
            } else {
              cat("***Need to refine bins for computing this shannon diversity of this GEDI metric***\n")
              cat(n_vals_out, "values outside the range of current bins \n")
              cat("The minimum is: ", min_val, "\n")
              cat("The maximum is: ", max_val, "\n") 
              cat("***\n")
            }

          } else {
            agg_stat <- rasterize(x = as.matrix(DT_r[, ..coords]), values = DT_r[,get(m)], y = base, fun = stats_dict[[s]])
          }
          names(agg_stat) <- paste0(m_lab, "_", as.character(names(stats_dict[s])))
          m_rast_list[[s]] <- agg_stat
          
          # check for valid values 
          m_ck <-ifelse(uniqueN(x = values(agg_stat), na.rm = TRUE) > 0, TRUE, FALSE)
          m_ck_list[[s]] <- m_ck
          
        } # ends the stats loop
        rm(DT_r)

        # check that at least 1 stat raster has valid values
        if (any(unlist(m_ck_list))){ 
          # stack all of the statistics for the metric of interest and export
          agg_stats <- rast(m_rast_list)
          shotctf_mask <- agg_stats[paste0(m_lab, "_countf")]
          shotctf_mask[shotctf_mask <= 1] <- NA
          
          # mask out pixels that have 1 shot or less
          agg_stats_masked = mask(agg_stats, shotctf_mask)
          
          cat(">> Saving", m_lab, "stats raster \n")
          writeRaster(x = crop(agg_stats_masked, ext_crop), filename = paste0(save_dir, prefix, m_lab, "_vf_", date_dec_start_lab, "_", date_dec_end_lab, "_", as.character(res), "m.tif"), 
                      gdal=c("COMPRESS=LZW", "TILED=YES"), NAflag = -9999, overwrite = TRUE)
          cat("\n")
          rm(agg_stat, m_rast_list, agg_stats, agg_stats_masked)
        } else {
          cat(">> NOT saving", m_lab, "stats raster because there are 0 valid values \n") 
          cat("\n")
        }
      }# ends the metric loop
      rm(agg_shotct, agg_shotct_g, shotct_mask, shotct_mask_g)
      
    } else { # if there are no shots within this time period
      cat("  there are 0 high quality veg. shots to rasterize during this time period \n")
    }
    rm(DT_tv)
    
    #Increase the date index counter
    cat(" ~~~ End Date Range: ", date_dec_start_lab, "to", date_dec_end_lab, "~~~\n")
    di <- di + 2
    cat("\n\n")
    
  } # ends the date pair loop
  cat("----- End rasterization at", res, "m spatial resolution ----- \n")
  cat("\n\n\n")
  
}# ends the base raster loop
cat("End of gedi_l2-l4a_s03_rasterize_chunks.R")
cat("Warnings from the R script: \n")
warnings()
cat("End of warnings \n")
cat("\n")

