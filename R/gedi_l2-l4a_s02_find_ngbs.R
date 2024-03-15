# Title: gedi_l2-l4a_s02_find_ngbs.R
# Author: Patrick Burns [pb463@nau.edu]
# Purpose: find neighbors of each chunk
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

# 1. directory containing summaries of GEDI processed orbits
summ_comb_tab <- as.character(args[1])

# 2. distance to search from each chunk
chunk_add_dist <- as.numeric(args[2])

# 3. directory to save files which list neighbors
save_dir <- as.character(args[3])



#####
# 2. Processing
#####
gedi_comb_chunk_summs(summ_comb_tab = summ_comb_tab,
                      chunk_add_dist = chunk_add_dist,
                      save_dir = save_dir)
