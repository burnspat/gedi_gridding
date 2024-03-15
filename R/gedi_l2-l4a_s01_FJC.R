# Title: gedi_l2-l4a_s01_FJC.R
# Author: Patrick Burns [pb463@nau.edu]
# Purpose: filter, join, and clip GEDI L2 and L4A data. To be run in conjunction with bash script which uses SLURM Arrays
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

# 1.2.3. L2A, L2B, L4A input file paths
l2a_input <- args[1]
l2b_input <- args[2]
l4a_input <- args[3]

# 4. name of the region
reg_name <- args[4]

# 5. region geometry
reg_geom <- args[5]

# 6. chunking dimension
chunk_dim <- as.numeric(args[6])

# 7. save directory
save_dir <- args[7]



#####
# 2. Processing
#####

cat("--- L2A data extraction and filtering --- \n")
cat("\n")
# L2A extraction
g2a_e <- gedi_l2a002_extract(l2a_input = l2a_input)
# L2A filtering
g2a_f <- gedi_l2a002_filter(l2_DT = g2a_e,
                            gf_algo = 0,
                            apply_l2a_qflag = TRUE,
                            apply_deg_flag = FALSE,
                            min_sens = 0.95,
                            min_sens_tegb = 0.98,
                            max_dem_diff = 150,
                            cover_night_only = FALSE,
                            leaf_on_only = FALSE,
                            max_ls_waterp = NULL,
                            max_urb_prop = NULL,
                            save_file = NULL)
rm(g2a_e)
cat("\n")


cat("--- L2B data extraction and filtering --- \n")
cat("\n")
# L2B extraction
g2b_e <- gedi_l2b002_extract(l2b_input = l2b_input)
# L2B filtering
# note: filter using low sensitivity since shots will be joined to L2A in the next step
g2b_f <- gedi_l2b002_filter(l2_DT = g2b_e, 
                            apply_algrun_flag = FALSE,
                            apply_l2a_qflag = TRUE,
                            apply_l2b_qflag = FALSE,                 
                            apply_deg_flag = FALSE,
                            max_dem_diff = 150,
                            cover_night_only = FALSE,
                            leaf_on_only = FALSE,
                            max_ls_waterp = NULL,
                            max_urb_prop = NULL,
                            save_file = NULL)
rm(g2b_e)
cat("\n")


cat("--- L2A and L2B joining --- \n")
cat("\n")
g2_j <- gedi_l2_join(l2a_DT = g2a_f, l2b_DT = g2b_f)
rm(g2a_f, g2b_f)
cat("\n")


cat("--- L4A data extraction and filtering --- \n")
# L4A extraction
g4a_e <- gedi_l4a002_extract(l4a_input = l4a_input)
# L4A filtering
# note: filter using low sensitivity since shots will be joined to L2 in the next step
g4a_f <- gedi_l4a002_filter(l4a_DT = g4a_e,
                            apply_algrun_flag = FALSE,
                            apply_l2a_qflag = TRUE,
                            apply_l4a_qflag = FALSE,
                            apply_deg_flag = FALSE,
                            min_sens = 0.95,
                            min_sens_tegb = 0.98,
                            cover_night_only = FALSE,
                            leaf_on_only = FALSE,
                            max_ls_waterp = NULL,
                            max_urb_prop = NULL,
                            save_file = NULL)
rm(g4a_e)
cat("\n")


cat("--- L2 and L4A joining --- \n")
g24a_j <- gedi_l2l4a_join(l2_DT = g2_j, l4a_DT = g4a_f)
rm(g2_j, g4a_f)
cat("\n")


# Exclude bad granules based on local outlier detection
# TODO: update bad gran list
cat("--- Identifying bad granules --- \n")
g24a_o <- gedi_outliers_umd(table = g24a_j,
                            grid_shp = '/projects/above_gedi/users/pburns/GEDI/grids/grid_EASE2_72km_land_pm52_named.shp',
                            grid_name = 'EASE72_id',
                            outliers_tab = '/projects/above_gedi/users/pburns/GEDI/bad_granules/issgedi_l4b_excluded_granules_r002_thruMW222_rec20230917.csv',
                            remove = FALSE,
                            save_file = NULL)
rm(g24a_j)
cat("\n")

# Add a high quality flag for vegetation and optionally remove low quality shots
cat("--- Adding a high quality flag --- \n")
g24a_q <- gedi_l2l4a_hqf(table = g24a_o,
                         remove = FALSE,
                         save_file = NULL)
rm(g24a_o)
cat("\n")

# compute custom metrics
cat("--- Computing custom metrics --- \n")
g24a_m <- gedi_l2l4a002_moremets(l24a_DT = g24a_q,
                                 save_file = NULL)
rm(g24a_q)
cat("\n")

# Chunk the GEDI shots into regular grids
cat("--- Chunking shots --- \n")
save_table <- paste0(save_dir, 'tables/',
                    tools::file_path_sans_ext(sub(basename(l2a_input), pattern = "GEDI02_A", replacement = "gedi0204a")),
                    "_filt.csv")
summ_save_table <- paste0(save_dir, 'summaries/',
                         tools::file_path_sans_ext(sub(basename(l2a_input), pattern = "GEDI02_A", replacement = "gedi0204a")),
                         "_summ.csv")
gedi_chunk_gcs(table = g24a_m, chunk_dim = chunk_dim, chunk_x_max = 180, chunk_y_max = 52.0, 
                save_file = save_table, summ_file = summ_save_table)