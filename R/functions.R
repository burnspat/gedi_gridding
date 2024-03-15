# Libraries
.libPaths(c("/projects/above_gedi/users/pburns/envs/geospatial/lib/R/library"))
library(sf, warn.conflicts = F, quietly = T)
library(data.table, warn.conflicts = F, quietly = T)
library(rhdf5, warn.conflicts = F, quietly = T)
library(raster, warn.conflicts = F, quietly = T)
library(terra, warn.conflicts = F, quietly = T)
library(sp, warn.conflicts = F, quietly = T)
library(lubridate, warn.conflicts = F, quietly = T)
library(httr, warn.conflicts = F, quietly = T)
library(spatstat, warn.conflicts = F, quietly = T)

# Test files
# l2a_input <- '/scratch/pb463/temp/gedi_l2a/GEDI02_A_2020193060840_O08943_02_T00094_02_003_01_V002.h5'
# l2b_input <- '/scratch/pb463/temp/gedi_l2b/GEDI02_B_2019171151005_O02947_02_T01823_02_003_01_V002.h5'
# borneo_shp <- '/scratch/pb463/temp/borneo.shp'



# Function: gedi_finder
# Author: Cole Krehbiel, LPDAAC
# Last Updated: 05/10/2021
# Purpose: spatially query NASA CMR for GEDI L1B, L2A, or L2B granules
# Arguments:
# > product = "GEDI02_A.002" or "GEDI02_B.002"
# > bbox = Bounding box corner coordinates - LL Longitude, LL Latitude, UR Longitude, UR Latitude
# Returns: list of granules
gedi_finder <- function(product = NULL,
                        bbox = NULL) {

  # Define the base CMR granule search url, including LPDAAC provider name and max page size (2000 is the max allowed)
  if (product == "GEDI01_B.002" | product == "GEDI02_A.002" | product == "GEDI02_B.002"){
    cmr <- "https://cmr.earthdata.nasa.gov/search/granules.json?pretty=true&provider=LPDAAC_ECS&page_size=2000&concept_id="
  } else if (product == "GEDI04_A.002"){
    cmr <- "https://cmr.earthdata.nasa.gov/search/granules.json?pretty=true&provider=ORNL_CLOUD&page_size=2000&concept_id="
  }

  # Set up dictionary where key is GEDI shortname + version and value is CMR Concept ID
  concept_ids <- list('GEDI01_B.002'='C1908344278-LPDAAC_ECS',
                      'GEDI02_A.002'='C1908348134-LPDAAC_ECS',
                      'GEDI02_B.002'='C1908350066-LPDAAC_ECS',
                      'GEDI04_A.002'='C2237824918-ORNL_CLOUD')


  # CMR uses pagination for queries with more features returned than the page size
  page <- 1
  bbox <- sub(' ', '', bbox)  # Remove any white spaces
  granules <- list()          # Set up a list to store and append granule links to

  # Send GET request to CMR granule search endpoint w/ product concept ID, bbox & page number
  cmr_response <- GET(sprintf("%s%s&bounding_box=%s&pageNum=%s", cmr, concept_ids[[product]],bbox,page))

  # Verify the request submission was successful
  if (cmr_response$status_code==200){

    # Send GET request to CMR granule search endpoint w/ product concept ID, bbox & page number, format return as a list
    cmr_response <- content(GET(sprintf("%s%s&bounding_box=%s&pageNum=%s", cmr, concept_ids[[product]],bbox,page)))$feed$entry

    # If 2000 features are returned, move to the next page and submit another request, and append to the response
    while(length(cmr_response) %% 2000 == 0){
      page <- page + 1
      cmr_response <- c(cmr_response, content(GET(sprintf("%s%s&bounding_box=%s&pageNum=%s", cmr, concept_ids[[product]],bbox,page)))$feed$entry)
    }

    # CMR returns more info than just the Data Pool links, below use for loop to go through each feature, grab DP link, and add to list
    for (i in 1:length(cmr_response)) {
      granules[[i]] <- cmr_response[[i]]$links[[1]]$href
    }

    # Return the list of links
    return(granules)
  } else {

    # If the request did not complete successfully, print out the response from CMR
    print(content(GET(sprintf("%s%s&bounding_box=%s&pageNum=%s", cmr, concept_ids[[product]],bbox,page)))$errors)
  }
} # End gedi_finder function



# Function: gedi_l2a002_extract
# Author: Patrick Burns, NAU
# Last Updated: 07/07/2023
# Purpose: extract and filter shots from GEDI L2A orbit files
# Arguments:
# > l2a_input = L2A .h5 path/file to process
# > gf_algos = vector of L2A ground finding algorithms to export. 0 for default. Other options: 1,2,3,4,5,6
# > rh_vals = vector of relative height percentiles to extract
# > save_file = where to save the output csv. If NULL, nothing will be saved
# Returns: data.table of filtered GEDI shots and selected RH profiles
gedi_l2a002_extract <- function(l2a_input = NULL,
                                gf_algos = c(0),
                                rh_vals = c(seq(0,90,5),92, 95, 98, 99, 100),
                                save_file = NULL){

  cat("Extracting observations from", l2a_input, "\n")
  cat("\n")

  # Define beams
  beams_pow = paste0("BEAM", c("0101", "0110", "1000", "1011"))
  beams_cov = paste0("BEAM", c("0000", "0001", "0010", "0011"))
  beams = c(beams_pow, beams_cov)

  # Create an empty list to store extracted observations
  datalist = list()

  # Loop through beams and extract shots
  t_s = Sys.time()
  for (b in 1:length(beams)){
    cat("Starting beam", b, "\n")

    # Read beam from file
    h = rhdf5::h5read(file = l2a_input, name = beams[b], bit64conversion='bit64')

    # Check to make sure the beam has observations
    totalObs = nrow(h$beam)

    if(!totalObs > 0){
      # Return an empty list if no observations
      datalist[[b]] = list()

      cat("Finished beam", b, "-", 0, "footprints \n")

    } else if(totalObs > 0){

      # There are some observations...
      # First work on time
      delta_time = h$delta_time
      date0 = lubridate::ymd_hms("2018-01-01 00:00:00")
      date_time = date0 + lubridate::seconds(delta_time)
      date = lubridate::date(date_time)
      date_decimal = lubridate::decimal_date(date_time)
      doy = lubridate::yday(date_time)
      year = lubridate::year(date_time)
      month = lubridate::month(date_time)
      day = lubridate::day(date_time)
      hour = lubridate::hour(date_time)
      min = lubridate::minute(date_time)
      sec = lubridate::second(date_time)
      time = substr(date_time,12,19)
      
      # Make a data.table to hold general shot data and landcover data
      DT = data.table(
        # general shot info
        shot_num = as.character(h$shot_number), # TODO: confirm UINT64 is ok
        orbit = as.character(unlist(strsplit(x = basename(l2a_input), split = "_"))[4]),
        sub_orbit = as.numeric(unlist(strsplit(x = basename(l2a_input), split = "_"))[5]),
        beam = as.character(rep(beams[b], nrow(h$shot_number))),
        delta_time = as.numeric(delta_time),
        date = as.character(date),
        date_dec = as.numeric(date_decimal),
        time = as.character(time),
        doy = as.integer(doy),
        sol_elev = as.numeric(h$solar_elevation),
        sol_azim = as.numeric(h$solar_azimuth),
        lon_lm_a0 = as.numeric(h$lon_lowestmode),
        lat_lm_a0 = as.numeric(h$lat_lowestmode),
        elev_lm_a0 = as.numeric(h$elev_lowestmode),
        elev_hr_a0 = as.numeric(h$elev_highestreturn),
        elev_dem_srtm = as.numeric(h$digital_elevation_model_srtm),
        elev_dem_tdx = as.numeric(h$digital_elevation_model),
        stale_flag = as.integer(h$geolocation$stale_return_flag),
        deg_flag = as.integer(h$degrade_flag),
        surf_flag = as.integer(h$surface_flag),
        sens_a0 = as.numeric(h$sensitivity),
        l2a_qflag_a0 = as.integer(h$quality_flag),
        l2a_qflag_a1 = as.integer(h$geolocation$quality_flag_a1),
        l2a_qflag_a2 = as.integer(h$geolocation$quality_flag_a2),
        l2a_qflag_a3 = as.integer(h$geolocation$quality_flag_a3),
        l2a_qflag_a4 = as.integer(h$geolocation$quality_flag_a4),
        l2a_qflag_a5 = as.integer(h$geolocation$quality_flag_a5),
        l2a_qflag_a6 = as.integer(h$geolocation$quality_flag_a6),
        #l2a_qflag_a10 = as.integer(h$geolocation$quality_flag_a10),
        l2a_selalg_a0 = as.integer(h$selected_algorithm),
        energy_total = as.numeric(h$energy_total),
        num_modes_a0 = as.integer(h$num_detectedmodes),

        # landcover
        ls_treecov = as.numeric(h$land_cover_data$landsat_treecover),
        ls_waterp = as.integer(h$land_cover_data$landsat_water_persistence),
        leafoff_flag = as.integer(h$land_cover_data$leaf_off_flag),
        leafoff_doy = as.integer(h$land_cover_data$leaf_off_doy),
        leafon_doy = as.integer(h$land_cover_data$leaf_on_doy),
        pft = as.integer(h$land_cover_data$pft_class),
        region = as.integer(h$land_cover_data$region_class),
        urb_prop = as.integer(h$land_cover_data$urban_proportion)
      )


      # Extract RH profiles associated with different ground-finding algos
      rh_gf_list <- list()
      gf_counter <- 1

      for (g in gf_algos){
        cat("Working on GF Algorithm ", g, "\n")
        if (g == 0){
          # This is the default ground finding algorithm
          rh_all <- t(h$rh)
          rh_sel <- as.data.table(rh_all[,rh_vals + 1])
          setnames(x = rh_sel, old = colnames(rh_sel), new = paste0("rh_",rh_vals,"_a",g))

          # If only using default algo, grab sensitivity values from other gf_algos
          if (length(g)==1 && g==0){
            rh_sel$sens_a1 <- h$geolocation$sensitivity_a1
            rh_sel$sens_a2 <- h$geolocation$sensitivity_a2
            rh_sel$sens_a3 <- h$geolocation$sensitivity_a3
            rh_sel$sens_a4 <- h$geolocation$sensitivity_a4
            rh_sel$sens_a5 <- h$geolocation$sensitivity_a5
            rh_sel$sens_a6 <- h$geolocation$sensitivity_a6
          }

          # Store the result in a list
          rh_gf_list[[gf_counter]] <- rh_sel

          gf_counter <- gf_counter + 1
        } else {
          keep_1d_base <- c("lon_lowestmode", "lat_lowestmode", "elev_lowestmode",
                            "sensitivity","quality_flag", "num_detectedmodes")
          change_arr <- c("lon_lowestmode", "lat_lowestmode", "elev_lowestmode",
                          "sensitivity")
          change_chr <- c("quality_flag", "num_detectedmodes")
          keep_1d <- paste0(keep_1d_base, "_a", g)
          keep_1d_chr <- paste0(change_chr, "_a", g)
          keep_1d_arr <- paste0(change_arr, "_a", g)
          geo_sel <- as.data.table(h$geolocation[keep_1d])

          #Change column types
          geo_sel <- geo_sel[,(keep_1d_arr) := lapply(.SD, as.numeric),
                             .SDcols = keep_1d_arr]
          geo_sel <- geo_sel[,(keep_1d_chr) := lapply(.SD, as.integer),
                             .SDcols = keep_1d_chr]
          setnames(x = geo_sel, old = colnames(geo_sel),
                   new = paste(c("lon_lm", "lat_lm", "elev_lm", "sens", "l2a_qflag", "num_modes"), g, sep = "_a"))

          # make a DT of the RH profile and convert to meters
          rh_get <- c(paste0("rh_a",g))
          rh_list <- h$geolocation[rh_get]
          rh_ext <- (t(rh_list[[1]])/100)[,rh_vals+1]
          rh_sel <- as.data.table(rh_ext)
          setnames(x = rh_sel, old = colnames(rh_sel), new = paste0("rh_",rh_vals,"_a",g))

          # Combine geolocation, quality, and RH profile for this GF algo
          rh_gf_list[[gf_counter]] <- cbind(geo_sel, rh_sel)

          gf_counter <- gf_counter + 1
        }
      }#end ground-finding algo loop

      # Merge RH profiles for all GF algos
      DT_RH <- do.call(cbind, rh_gf_list)
      cat("Combined RH profiles associated with selected ground-finding algorithms \n")
    }

    # Merge inital DT and DT_RH
    DT_comb <- cbind(DT, DT_RH)
    cat("Combined general shot info with RH profiles \n")

    # Store the extracted and filtered shots in a list
    datalist[[b]] <- DT_comb

    h5closeAll()

    cat("Finished beam", b, "-", nrow(DT_comb), "footprints \n")
    cat("\n")

  } # End beam loop extraction

  cat("Extracted info from L2A HDF5 files \n")

  # Merge list into one data.table
  combined = do.call(rbind, datalist)
  cat("Found a combined total of", nrow(combined), "footprints \n")
  cat("\n")

  # Check total time to read hdf and make df
  t_e = Sys.time()
  t_d = t_e - t_s
  cat("Time to process h5 file:", t_d, "\n")

  if (nrow(combined) == 0){
    cat("No shots found \n")
    cat("\n")
    quit(save = "no")
  } else {
    if (!is.null(save_file)){
      f_base <- tools::file_path_sans_ext(basename(l2a_input))
      fwrite(x = combined, row.names = FALSE, file = save_file)
      cat("Saved extracted footprints to", save_file, "\n")
      cat("\n")
    } else {
      # Return a data.table
      cat("Did NOT save extracted footprints to .csv \n")
      cat("\n")
      return(combined)
    }
  }
} # End gedi_l2a002_extract function



# Function: gedi_l2b002_extract
# Author: Patrick Burns, NAU
# Last Updated: 06/12/2023
# Purpose: extract and filter shots from GEDI L2B orbit files
# Arguments:
# > l2b_input = L2B .h5 path/file to process
# > vprofs = whether or not to export the vertical cover, PAI, and PAVD vertical profiles. TRUE/FALSE
# > save_file = where to save the output csv. If NULL, nothing will be saved
# Returns: data.table of filtered GEDI shots and selected vertical profile metrics
gedi_l2b002_extract <- function(l2b_input = NULL,
                                vprofs = TRUE,
                                save_file = NULL){

  cat("Extracting observations from", l2b_input, "\n")
  cat("\n")

  # Define beams
  beams_pow = paste0("BEAM", c("0101", "0110", "1000", "1011"))
  beams_cov = paste0("BEAM", c("0000", "0001", "0010", "0011"))
  beams = c(beams_pow, beams_cov)

  # Create an empty list to store extracted observations
  datalist = list()

  # Loop through beams and extract shots
  t_s = Sys.time()
  for (b in 1:length(beams)){
    cat("Starting beam", b, "\n")

    # Read beam from file
    h = rhdf5::h5read(file = l2b_input, name = beams[b], bit64conversion='bit64')

    # Check to make sure the beam has observations
    totalObs = nrow(h$beam)

    if(!totalObs > 0){
      # Return an empty list if no observations
      datalist[[b]] = list()

      cat("Finished beam", b, "-", 0, "footprints \n")

    } else if(totalObs > 0){

      # There are some observations...
      # First work on time
      delta_time = h$geolocation$delta_time
      date0 = lubridate::ymd_hms("2018-01-01 00:00:00")
      date_time = date0 + lubridate::seconds(delta_time)
      date = lubridate::date(date_time)
      date_decimal = lubridate::decimal_date(date_time)
      doy = lubridate::yday(date_time)
      year = lubridate::year(date_time)
      month = lubridate::month(date_time)
      day = lubridate::day(date_time)
      hour = lubridate::hour(date_time)
      min = lubridate::minute(date_time)
      sec = lubridate::second(date_time)
      time = substr(date_time,12,19)

      # Make a data.table to hold general shot data and landcover data
      DT = data.table(
        # general shot info
        shot_num = as.character(h$shot_number), # TODO: confirm UINT64 is ok
        orbit = as.character(unlist(strsplit(x = basename(l2b_input), split = "_"))[4]),
        sub_orbit = as.numeric(unlist(strsplit(x = basename(l2b_input), split = "_"))[5]),
        beam = as.character(rep(beams[b], nrow(h$shot_number))),
        delta_time = as.numeric(delta_time),
        date = as.character(date),
        date_dec = as.numeric(date_decimal),
        time = as.character(time),
        doy = as.integer(doy),
        sol_elev = as.numeric(h$geolocation$solar_elevation),
        sol_azim = as.numeric(h$geolocation$solar_azimuth),
        lon_lm_a0 = as.numeric(h$geolocation$lon_lowestmode),
        lat_lm_a0 = as.numeric(h$geolocation$lat_lowestmode),
        elev_lm_a0 = as.numeric(h$geolocation$elev_lowestmode),
        elev_hr_a0 = as.numeric(h$geolocation$elev_highestreturn),
        elev_dem_tdx = as.numeric(h$geolocation$digital_elevation_model),
        stale_flag = as.integer(h$stale_return_flag), 
        deg_flag = as.integer(h$geolocation$degrade_flag),
        surf_flag = as.integer(h$surface_flag),
        beam_azim = as.numeric(h$geolocation$local_beam_azimuth),
        beam_elev = as.numeric(h$geolocation$local_beam_elevation),
        sens_a0 = as.numeric(h$sensitivity),
        l2a_qflag_a0 = as.integer(h$l2a_quality_flag),
        l2b_qflag_a0 = as.integer(h$l2b_quality_flag),
        l2b_algrun_flag = as.integer(h$algorithmrun_flag),
        l2a_sel_alg_a0 = as.integer(h$selected_l2a_algorithm),
        num_modes_a0 = as.integer(h$num_detectedmodes),
        rh_100_a0 = as.numeric(h$rh100)/100,
        cover_a0 = as.numeric(h$cover),
        pai_a0 = as.numeric(h$pai),
        fhd_pai_1m_a0 = as.numeric(h$fhd_normal),
        omega = as.numeric(h$omega),
        rossg = as.numeric(h$rossg),
        pgap_theta_a0 = as.numeric(h$pgap_theta),
        pgap_theta_err_a0 = as.numeric(h$pgap_theta_error),
        

        # landcover
        ls_treecov = as.numeric(h$land_cover_data$landsat_treecover),
        ls_waterp = as.integer(h$land_cover_data$landsat_water_persistence),
        leafoff_flag = as.integer(h$land_cover_data$leaf_off_flag),
        leafoff_doy = as.integer(h$land_cover_data$leaf_off_doy),
        leafon_doy = as.integer(h$land_cover_data$leaf_on_doy),
        pft = as.integer(h$land_cover_data$pft_class),
        region = as.integer(h$land_cover_data$region_class),
        urb_prop = as.integer(h$land_cover_data$urban_proportion)
      )


      if (vprofs){
        # Extract vertical profiles

        # Function to rearrange vertical profile values so that cumulative profile starts from 0
        # rearr <- function(x){
        #   valid_ncol <- tail(x = x, n = 1)
        #   if (is.na(valid_ncol)){
        #     flip <- rep(-9999,30)
        #   } else if(valid_ncol == 0) {
        #     flip <- rep(-9999,30)
        #   } else {
        #     flip <- c(rev(x[1:valid_ncol]), rep(0,30-valid_ncol))
        #   }
        #   return(flip)
        # }

        # cover profile extraction
        prof_sel_i_cov <- as.data.table(t(h[["cover_z"]]))

        # cover_z is cumulative (high to low) starting from column 1
        #prof_sel_i_cov <- prof_sel_i_cov[, valid_ncol := rowSums(prof_sel_i_cov > 0)]

        # Rearrange values in each row so they start from 0
        # prof_sel_r_cov <- as.data.table(t(apply(X = prof_sel_i_cov, MARGIN = 1, FUN = rearr)))

        # Set names
        prof_new_names = paste0('cover_l', seq(1,30))
        setnames(x = prof_sel_i_cov, old = colnames(prof_sel_i_cov[,1:30]), new = prof_new_names)

        # Only keep cover up to 100 m
        prof_sel_cov <- prof_sel_i_cov[,1:21]
        rm(prof_sel_i_cov)


        # PAI profile extraction
        prof_sel_i_pai <- as.data.table(t(h[["pai_z"]]))

        # pai_z is cumulative (high to low) starting from column 1
        # prof_sel_i_pai <- prof_sel_i_pai[, valid_ncol := rowSums(prof_sel_i_pai > 0)]

        # Rearrange values in each row so they start from 0
        # prof_sel_r_pai <- as.data.table(t(apply(X = prof_sel_i_pai, MARGIN = 1, FUN = rearr)))

        # Set names
        prof_new_names = paste0('pai_l', seq(1,30))
        setnames(x = prof_sel_i_pai, old = colnames(prof_sel_i_pai[,1:30]), new = prof_new_names)

        # Only keep pai up to 100 m
        prof_sel_pai <- prof_sel_i_pai[,1:21]
        rm(prof_sel_i_pai)


        # PAVD profile extraction
        prof_sel_pav <- as.data.table(t(h[["pavd_z"]]))
        setnames(x = prof_sel_pav, old = colnames(prof_sel_pav), new = paste0("pavd_", seq(0,5*29,5), "_", seq(5,5*30,5)))

        # Only keep pavd up to 100 m
        prof_sel_pav <- prof_sel_pav[,1:20]

        # Merge vertical profile metrics
        DT_VP <- cbind(prof_sel_cov, prof_sel_pai, prof_sel_pav)
        cat("Combined vertical profile metrics \n")

        # Merge inital DT and DT_VP
        DT_comb <- cbind(DT, DT_VP)
        cat("Combined general shot info with vertical profile metrics \n")
      } else {
        # Don't extract vertical profiles
        DT_comb <- DT
      }

      # Store the extracted and filtered shots in a list
      datalist[[b]] <- DT_comb

      h5closeAll()

      cat("Finished beam", b, "-", nrow(DT_comb), "footprints \n")
      cat("\n")
    }

  } # End beam loop extraction

  cat("Extracted info from L2B HDF5 files \n")

  # Merge list into one data.table
  combined = do.call(rbind, datalist)
  cat("Found a combined total of", nrow(combined), "footprints \n")
  cat("\n")

  # Check total time to read hdf and make df
  t_e = Sys.time()
  t_d = t_e - t_s
  cat("Time to process h5 file:", t_d, "\n")

  if (nrow(combined) == 0){
    cat("No shots found \n")
    cat("\n")
    quit(save = "no")
  } else {
    if (!is.null(save_file)){
      f_base <- tools::file_path_sans_ext(basename(l2b_input))
      fwrite(x = combined, row.names = FALSE, file = save_file)
      cat("Saved extracted footprints to", save_file, "\n")
      cat("\n")
    } else {
      # Return a data.table
      cat("Did NOT save extracted footprints to .csv \n")
      cat("\n")
      return(combined)
    }
  }
} # End gedi_l2b002_extract function



# Function: gedi_l4a002_extract
# Author: Patrick Burns, NAU
# Last Updated: 06/12/2023
# Purpose: extract and filter shots from GEDI L4A orbit files
# Arguments:
# > l4a_input = L4A .h5 path/file to process
# > save_file = where to save the output csv. If NULL, nothing will be saved
# Returns: data.table of filtered GEDI L4A shots
gedi_l4a002_extract <- function(l4a_input = NULL,
                                save_file = NULL){

  cat("Extracting observations from", l4a_input, "\n")
  cat("\n")

  # Define beams
  beams_pow = paste0("BEAM", c("0101", "0110", "1000", "1011"))
  beams_cov = paste0("BEAM", c("0000", "0001", "0010", "0011"))
  beams = c(beams_pow, beams_cov)

  # Create an empty list to store extracted observations
  datalist = list()

  # Loop through beams and extract shots
  t_s = Sys.time()
  for (b in 1:length(beams)){
    cat("Starting beam", b, "\n")

    # Read beam from file
    h = rhdf5::h5read(file = l4a_input, name = beams[b], bit64conversion='bit64')

    # Check to make sure the beam has observations
    totalObs = nrow(h$beam)

    if(!totalObs > 0){
      # Return an empty list if no observations
      datalist[[b]] = list()

      cat("Finished beam", b, "-", 0, "(", 0, "% )", "footprints \n")

    } else if(totalObs > 0){

      # There are some observations...
      # First work on time
      delta_time = h$delta_time
      date0 = lubridate::ymd_hms("2018-01-01 00:00:00")
      date_time = date0 + lubridate::seconds(delta_time)
      date = lubridate::date(date_time)
      date_decimal = lubridate::decimal_date(date_time)
      doy = lubridate::yday(date_time)
      year = lubridate::year(date_time)
      month = lubridate::month(date_time)
      day = lubridate::day(date_time)
      hour = lubridate::hour(date_time)
      min = lubridate::minute(date_time)
      sec = lubridate::second(date_time)
      time = substr(date_time,12,19)

      # Make a data.table to hold general shot data and landcover data
      DT = data.table(
        # general shot info
        shot_num = as.character(h$shot_number), # TODO: confirm UINT64 is ok
        orbit = as.character(unlist(strsplit(x = basename(l4a_input), split = "_"))[4]),
        sub_orbit = as.numeric(unlist(strsplit(x = basename(l4a_input), split = "_"))[5]),
        beam = as.character(rep(beams[b], nrow(h$shot_number))),
        delta_time = as.numeric(delta_time),
        date = as.character(date),
        date_dec = as.numeric(date_decimal),
        time = as.character(time),
        doy = as.integer(doy),
        sol_elev = as.numeric(h$solar_elevation),
        lon_lm_a0 = as.numeric(h$lon_lowestmode),
        lat_lm_a0 = as.numeric(h$lat_lowestmode),
        elev_lm_a0 = as.numeric(h$elev_lowestmode),
        stale_flag = as.integer(h$geolocation$stale_return_flag),
        deg_flag = as.integer(h$degrade_flag),
        surf_flag = as.integer(h$surface_flag),
        sens_a0 = as.numeric(h$sensitivity),
        sens_a1 = as.numeric(h$geolocation$sensitivity_a1),
        sens_a2 = as.numeric(h$geolocation$sensitivity_a2),
        sens_a3 = as.numeric(h$geolocation$sensitivity_a3),
        sens_a4 = as.numeric(h$geolocation$sensitivity_a4),
        sens_a5 = as.numeric(h$geolocation$sensitivity_a5),
        sens_a6 = as.numeric(h$geolocation$sensitivity_a6),
        sens_a10 = as.numeric(h$geolocation$sensitivity_a10),
        l2a_qflag_a0 = as.integer(h$l2_quality_flag),
        l4a_qflag_a0 = as.integer(h$l4_quality_flag),
        l2a_selalg_a0 = as.integer(h$selected_algorithm),
        agbd_a0 = as.numeric(h$agbd),
        agbd_pi_low_a0 = as.numeric(h$agbd_pi_lower),
        agbd_pi_upper_a0 = as.numeric(h$agbd_pi_upper),
        agbd_se_a0 = as.numeric(h$agbd_se),
        l4a_algrun_flag = as.integer(h$algorithm_run_flag),
        l4a_pred_lim_flag = as.integer(h$predictor_limit_flag),
        l4a_resp_lim_flag = as.integer(h$response_limit_flag),
        l4a_pred_strat = as.character(h$predict_stratum),


        # landcover
        ls_treecov = as.numeric(h$land_cover_data$landsat_treecover),
        ls_waterp = as.integer(h$land_cover_data$landsat_water_persistence),
        leafoff_flag = as.integer(h$land_cover_data$leaf_off_flag),
        leafoff_doy = as.integer(h$land_cover_data$leaf_off_doy),
        leafon_doy = as.integer(h$land_cover_data$leaf_on_doy),
        pft = as.integer(h$land_cover_data$pft_class),
        region = as.integer(h$land_cover_data$region_class),
        urb_prop = as.integer(h$land_cover_data$urban_proportion)
      )

      cat("Combined general shot info with AGBD metrics \n")  

      # Store the extracted and filtered shots in a list
      datalist[[b]] <- DT

      h5closeAll()

      cat("Finished beam", b, "-", nrow(DT), "footprints \n")
      cat("\n")
    }
  } # End beam loop extraction

  cat("Extracted info from L4A HDF5 files \n")

  # Merge list into one data.table
  combined = do.call(rbind, datalist)
  cat("Found a combined total of", nrow(combined), "footprints \n")
  cat("\n")

  # Check total time to read hdf and make df
  t_e = Sys.time()
  t_d = t_e - t_s
  cat("Time to process h5 file:", t_d, "\n")

  if (nrow(combined) == 0){
    cat("No shots found \n")
    cat("\n")
    quit(save = "no")
  } else {
    if (!is.null(save_file)){
      f_base <- tools::file_path_sans_ext(basename(l4a_input))
      fwrite(x = combined, row.names = FALSE, file = save_file)
      cat("Saved extracted footprints to", save_file, "\n")
      cat("\n")
    } else {
      # Return a data.table
      cat("Did NOT save extracted footprints to .csv \n")
      cat("\n")
      return(combined)
    }
  }
} # End gedi_l4a002_extract function



# Function: gedi_l2a002_filter
# Author: Patrick Burns, NAU
# Last Updated: 06/12/2023
# Purpose: filter GEDI L2A shots using various quality criteria
# Arguments:
# > l2_DT = extracted L2A/B shots as either a data.table or .csv file
# > gf_algo = ground-finding algorithm to use. Options: 0(default), 1,2,3,4,5
# > apply_l2a_qflag = whether or not to apply the L2A quality flag. TRUE/FALSE
# > apply_deg_flag - whether or not to apply the degrade flag filter. TRUE/FALSE
# > min_sens = the minimum sensitivity threshold to use for non-tropical evergreen broadleaf forests
# > min_sens_tegb = the minimum sensitivity threshold to use for tropical evergreen broadleaf forests
# > max_dem_diff - maximum allowable absolute elevation difference from TanDEM-X DEM
# > cover_night_only = whether or not to save coverage beam data acquired during the day
# > leaf_on_only = whether or not to only use leaf on shots
# > max_ls_waterp = maximum landsat water persistence value
# > max_urb_prop = maximum urban proportion value
# > save_file = where to save the output csv. If NULL, nothing will be saved and a data.table will be returned
# Returns: data.table or saved .csv of filtered GEDI L2A shots
gedi_l2a002_filter <- function(l2_DT = NULL,
                               gf_algo = 0,
                               apply_l2a_qflag = TRUE,
                               apply_deg_flag = TRUE,
                               min_sens = 0.95,
                               min_sens_tegb = 0.98,
                               max_dem_diff = 150,
                               cover_night_only = FALSE,
                               leaf_on_only = FALSE,
                               max_ls_waterp = 10,
                               max_urb_prop = 50,
                               save_file = NULL){

  # Handle l2_DT input
  type = class(l2_DT)
  if ('data.table' %in% type){ # the data.table already exists
    DT <- l2_DT
  } else if ('character' %in% type){ # read in .csv as data.table
    DT <- fread(file = l2_DT)
  }

  # Specify coverage and power beam names
  beams_pow = paste0("BEAM", c("0101", "0110", "1000", "1011"))
  beams_cov = paste0("BEAM", c("0000", "0001", "0010", "0011"))

  # Select the RH profile associated with specified ground-finding algorithm
  # rh_cols = grep(pattern = paste0("rh_.*a", gf_algo), x = names(DT), value=TRUE)
  # sens_cols = grep(pattern = "sens_a", x = names(DT), value=TRUE)
  # DT_rh = cbind(DT[,1:23], DT[, ..sens_cols], DT[,25:32], elev_diff_dem=DT$elev_diff_dem, DT[, ..rh_cols])

  # Remove duplicate columns
  # DT_u <- DT[, .SD, .SDcols = unique(names(DT))]

  # Preliminary common filters
  ground_elev = paste0("elev_lm_a",gf_algo)
  rh100 = paste0("rh_100_a", gf_algo)
  DT_f = DT[sens_a0 > 0.9 & sens_a0 <= 1
            & sens_a2 > 0.9 & sens_a2 <= 1
            & surf_flag == 1
            & stale_flag == 0
            & get(rh100) >= 0 & get(rh100) < 120
            & get(ground_elev) > -200
            & get(ground_elev) < 9000
            ]


  # Optional advanced filters
  # Whether or not to apply the L2A quality flag
  if (apply_l2a_qflag){
    DT_f <- DT_f[l2a_qflag_a0 == 1]
  }

  # Apply separate sensitivity filters for tropical evergreen broadleaf forests in specific regions
  DT_f_nontrop = DT_f[!(pft == 2 & region %in% c(4,5,6))
                     & sens_a2 > min_sens & sens_a2 <= 1]

  DT_f_trop = DT_f[(pft == 2 & region %in% c(4,5,6))
                   & sens_a2 > min_sens_tegb & sens_a2 <= 1]

  DT_f = rbind(DT_f_nontrop, DT_f_trop)

  # apply the maximum allowable DEM difference
  # Compare TanDEM-X elevations and GEDI ground elevation
  DT_f = DT_f[, elev_diff_dem:=(elev_dem_tdx - get(ground_elev))]
  if (!is.null(max_dem_diff)){
    DT_f <- DT_f[elev_diff_dem < max_dem_diff & elev_diff_dem > (-1*max_dem_diff)]
  }

  # filter by max surface water percentage
  if (!is.null(max_ls_waterp)){
    DT_f <- DT_f[ls_waterp < max_ls_waterp]
  }

  # filter by max urban proportion
  if (!is.null(max_urb_prop)){
    DT_f <- DT_f[urb_prop < max_urb_prop]
  }

  # Apply the degrade flag
  if (apply_deg_flag){
    # Only accept these degrade flag values
    deg_ok <- c(0,3,8,10,13,18,20,23,28,30,33,38,40,43,48,60,63,68)
    DT_f <- DT_f[deg_flag %in% deg_ok]
  }

  # Apply the leaf-on filter
  if (leaf_on_only){
    DT_f = DT_f[!(leafoff_flag == 1)]
  }

  # Apply filter to only use coverage beams at night
  if (cover_night_only){
    DT_f_pow = DT_f[beam %in% beams_pow]
    DT_f_cov = DT_f[sol_elev < 0 & beam %in% beams_cov]
    DT_f = rbind(DT_f_pow, DT_f_cov)
  } 

  # Check to see if there are shots after filtering
  if (nrow(DT_f) == 0){
    cat("0 shots remain after filtering \n")
    cat("\n")
    quit(save = "no")
  } else {
    cat(nrow(DT_f), "L2A shots remain after filtering \n")
    cat("\n")

    # Whether or not to save a .csv
    if (!is.null(save_file)){
      # Save as csv
      fwrite(x = DT_f, row.names = FALSE, file = save_file)
      cat("Saved filtered shots to", save_file, "\n")
      cat("\n")
    } else {
      # Return a data.table
      cat("Saved filtered shots as data.table \n")
      cat("\n")
      return(DT_f)
    }
  }
} # End gedi_l2a002_filter function



# Function: gedi_l2b002_filter
# Author: Patrick Burns, NAU
# Last Updated: 09/22/2022
# Purpose: filter GEDI L2B shots using various quality criteria
# Arguments:
# > l2_DT = extracted L2A/B shots as either a data.table or .csv file
# > apply_algrun_flag = whether or not to apply the algorithm run flag. TRUE/FALSE
# > apply_l2a_qflag = whether or not to apply the L2A quality flag. TRUE/FALSE
# > apply_l2b_qflag = whether or not to apply the L2B quality flag. TRUE/FALSE
# > apply_deg_flag - whether or not to apply the degrade flag filter. TRUE/FALSE
# > max_dem_diff - maximum allowable absolute elevation difference from TanDEM-X DEM
# > cover_night_only = whether or not to save coverage beam data acquired during the day
# > leaf_on_only = whether or not to only use leaf on shots
# > max_ls_waterp = maximum landsat water persistence value
# > max_urb_prop = maximum urban proportion value
# > save_file = where to save the output csv. If NULL, nothing will be saved and a data.table will be returned
# Returns: data.table of filtered GEDI L2A shots
gedi_l2b002_filter <- function(l2_DT = NULL,
                               apply_algrun_flag = TRUE, 
                               apply_l2a_qflag = TRUE,
                               apply_l2b_qflag = TRUE,                 
                               apply_deg_flag = TRUE,
                               max_dem_diff = 150,
                               cover_night_only = FALSE,
                               leaf_on_only = FALSE,
                               max_ls_waterp = 10,
                               max_urb_prop = 50,
                               save_file = NULL){

  # Handle l2_DT input
  type = class(l2_DT)
  if ('data.table' %in% type){ # the data.table already exists
    DT <- l2_DT
  } else if ('character' %in% type){ # read in .csv as data.table
    DT <- fread(file = l2_DT)
  }

  # Specify coverage and power beam names
  beams_pow = paste0("BEAM", c("0101", "0110", "1000", "1011"))
  beams_cov = paste0("BEAM", c("0000", "0001", "0010", "0011"))

  # Preliminary common filters
  DT_f = DT[sens_a0 > 0.9 & sens_a0 <= 1
            & surf_flag == 1
            & stale_flag == 0
            & rh_100_a0 >= 0 & rh_100_a0 < 120
            & elev_lm_a0 > -200
            & elev_lm_a0 < 9000
            ]


  # Optional advanced filters
  # apply the algorithm run flag
  if(apply_algrun_flag){
    DT_f = DT_f[l2b_algrun_flag == 1]
  }

  # Whether or not to apply the L2A quality flag
  if (apply_l2a_qflag){
    DT_f <- DT_f[l2a_qflag_a0 == 1]
  }

  # Whether or not to apply the L2B quality flag
  if (apply_l2b_qflag){
    DT_f <- DT_f[l2b_qflag_a0 == 1]
  }

  # # Apply separate sensitivity filters for tropical evergreen broadleaf forests
  # DT_f_nontrop = DT_f[!(pft == 2 & region %in% c(4,5,6))
  #                     & sens_a0 >= min_sens & sens_a0 <= 1]

  # DT_f_trop = DT_f[(pft == 2 & region %in% c(4,5,6))
  #                  & sens_a0 >= min_sens_tegb & sens_a0 <= 1]

  # DT_f = rbind(DT_f_nontrop, DT_f_trop)

  # apply the maximum allowable DEM difference
  # Compare TanDEM-X elevations and GEDI ground elevation
  DT_f = DT_f[, elev_diff_dem:=(elev_dem_tdx - elev_lm_a0)]
  if (!is.null(max_dem_diff)){
    DT_f <- DT_f[elev_diff_dem < max_dem_diff & elev_diff_dem > (-1*max_dem_diff)]
  }

  # filter by max surface water percentage
  if (!is.null(max_ls_waterp)){
    DT_f <- DT_f[ls_waterp < max_ls_waterp]
  }

  # filter by max urban proportion
  if (!is.null(max_urb_prop)){
    DT_f <- DT_f[urb_prop < max_urb_prop]
  }

  # Apply the degrade flag
  if (apply_deg_flag){
    # Only accept these degrade flag values
    deg_ok <- c(0,3,8,10,13,18,20,23,28,30,33,38,40,43,48,60,63,68)
    DT_f <- DT_f[deg_flag %in% deg_ok]
  }

  # Apply the leaf-on filter
  if (leaf_on_only){
    DT_f = DT_f[!(leafoff_flag == 1)]
  }

  # Apply filter to only use coverage beams at night
  if (cover_night_only){
    DT_f_pow = DT_f[beam %in% beams_pow]
    DT_f_cov = DT_f[sol_elev < 0 & beam %in% beams_cov]
    DT_f = rbind(DT_f_pow, DT_f_cov)
  } 

  # Check to see if there are shots after filtering
  if (nrow(DT_f) == 0){
    cat("0 shots remain after filtering \n")
    cat("\n")
    quit(save = "no")
  } else {
    cat(nrow(DT_f), "L2B shots remain after filtering \n")
    cat("\n")

    # Whether or not to save a .csv
    if (!is.null(save_file)){
      # Save as csv
      fwrite(x = DT_f, row.names = FALSE, file = save_file)
      cat("Saved filtered shots to", save_file, "\n")
      cat("\n")
    } else {
      # Return a data.table
      cat("Saved filtered shots as data.table \n")
      cat("\n")
      return(DT_f)
    }
  }
} # End gedi_l2b002_filter function



# Function: gedi_l4a002_filter
# Author: Patrick Burns, NAU
# Last Updated: 09/22/2022
# Purpose: filter GEDI L4A shots using various quality criteria
# Arguments:
# > l4a_DT = extracted L4A shots as either a data.table or .csv file
# > apply_algrun_flag = whether or not to apply the algorithm run flag
# > apply_l2a_qflag = whether or not to apply the L2A quality flag. TRUE/FALSE
# > apply_l4a_qflag = whether or not to apply the L4A quality flag. TRUE/FALSE
# > apply_deg_flag - whether or not to apply the degrade flag filter
# > min_sens = the minimum sensitivity threshold to use
# > min_sens_tegb = the minimum sensitivity threshold to use for tropical evergreen broadleaf forests
# > cover_night_only = whether or not to save coverage beam data acquired during the day
# > leaf_on_only = whether or not to only use leaf on shots
# > max_ls_waterp = maximum landsat water persistence value
# > max_urb_prop = maximum urban proportion value
# > save_file = where to save the output csv. If NULL, nothing will be saved and a data.table will be returned
# Returns: data.table of filtered GEDI L4A shots
gedi_l4a002_filter <- function(l4a_DT = NULL,
                               apply_algrun_flag = FALSE,
                               apply_l2a_qflag = TRUE,
                               apply_l4a_qflag = TRUE,
                               apply_deg_flag = FALSE,
                               min_sens = 0.95,
                               min_sens_tegb = 0.98,
                               cover_night_only = FALSE,
                               leaf_on_only = FALSE,
                               max_ls_waterp = NULL,
                               max_urb_prop = NULL,
                               save_file = NULL){

  # Handle l4a_DT input
  type = class(l4a_DT)
  if ('data.table' %in% type){ # the data.table already exists
    DT <- l4a_DT
  } else if ('character' %in% type){ # read in .csv as data.table
    DT <- fread(file = l4a_DT)
  }

  # Specify coverage and power beam names
  beams_pow <- paste0("BEAM", c("0101", "0110", "1000", "1011"))
  beams_cov <- paste0("BEAM", c("0000", "0001", "0010", "0011"))

  # Preliminary common filters
  DT_f = DT[sens_a0 > 0.9 & sens_a0 <= 1
            & sens_a2 > 0.9 & sens_a2 <= 1
            & surf_flag == 1
            & stale_flag == 0
            & elev_lm_a0 > -200
            & elev_lm_a0 < 9000
            ]
  
  # Optional advanced filters
  # apply the algorithm run flag
  if(apply_algrun_flag){
    DT_f = DT_f[l4a_algrun_flag == 1]
  }

  # Whether or not to apply the L2A quality flag
  if (apply_l2a_qflag){
    DT_f <- DT_f[l2a_qflag_a0 == 1]
  }

  # Whether or not to apply the L4A quality flag
  if (apply_l4a_qflag){
    DT_f <- DT_f[l4a_qflag_a0 == 1]
  }

  # Apply separate sensitivity filters for tropical evergreen broadleaf forests
  DT_f_nontrop = DT_f[!(pft == 2 & region %in% c(4,5,6))
                      & sens_a2 >= min_sens & sens_a2 <= 1]

  DT_f_trop = DT[(pft == 2 & region %in% c(4,5,6))
                 & sens_a2 >= min_sens_tegb & sens_a2 <= 1]

  DT_sens_f = rbind(DT_f_nontrop, DT_f_trop)
  

  # filter by max surface water percentage
  if (!is.null(max_ls_waterp)){
    DT_f <- DT_f[ls_waterp < max_ls_waterp]
  }

  # filter by max urban proportion
  if (!is.null(max_urb_prop)){
    DT_f <- DT_f[urb_prop < max_urb_prop]
  }
  
  # Apply the degrade flag
  if (apply_deg_flag){
    # Only accept these degrade flag values
    deg_ok <- c(0,3,8,10,13,18,20,23,28,30,33,38,40,43,48,60,63,68)
    DT_f <- DT_f[deg_flag %in% deg_ok]
  }

  # Apply the leaf-on filter
  if (leaf_on_only){
    DT_f = DT_f[!(leafoff_flag == 1)]
  }

  # Apply filter to only use coverage beams at night
  if (cover_night_only){
    DT_f_pow = DT_f[beam %in% beams_pow]
    DT_f_cov = DT_f[sol_elev < 0 & beam %in% beams_cov]
    DT_f = rbind(DT_f_pow, DT_f_cov)
  } 


  # Check to see if there are shots after filtering
  if (nrow(DT_f) == 0){
    cat("0 L4A shots remain after filtering \n")
    cat("\n")
    quit(save = "no")
  } else {
    cat(nrow(DT_f), "L4A shots remain after filtering \n")
    cat("\n")

    # Whether or not to save a .csv
    if (!is.null(save_file)){
      # Save as csv
      fwrite(x = DT_f, row.names = FALSE, file = save_file)
      cat("Saved filtered shots to", save_file, "\n")
      cat("\n")
    } else {
      # Return a data.table
      cat("Saved filtered shots as data.table \n")
      cat("\n")
      return(DT_f)
    }
  }
} # End gedi_l4a002_filter function



# Function: gedi_l24a002_moremets
# Author: Patrick Burns, NAU
# Last Updated: 06/29/2023
# Purpose: compute additional derived metrics from L2 + L4a
# Arguments:
# > l24a_DT = extracted L24A shots as either a data.table or .csv file
# > save_file = where to save the output csv. If NULL, nothing will be saved and a data.table will be returned
# Returns: data.table of GEDI L24A shots with additional metrics added
gedi_l2l4a002_moremets <- function(l24a_DT = NULL,
                                   save_file = NULL){

  # Handle l24a_DT input
  type = class(l24a_DT)
  if ('data.table' %in% type){ # the data.table already exists
    DT <- l24a_DT
  } else if ('character' %in% type){ # read in .csv as data.table
    DT <- fread(file = l24a_DT)
  }

  # get the original PAVD column names
  pavd_cols <- names(DT)[names(DT) %like% "pavd_"]

  # fraction of PAVD in the lowest height bin - pavd_0_5_frac
  DT_m <- DT[, pavd_0_5_frac:=fifelse(rh_100_a0 > 5, pavd_0_5/rowSums(.SD), 1), .SDcols = pavd_cols]

  # height of maximum PAVD
  DT_m <- DT_m[, pavd_max_h:=(max.col(.SD)*5), .SDcols = pavd_cols]

  # verictal distribution ratios of the RH profile - rhvdr_b, rhvdr_m, rhvdr_t
  DT_m <- DT_m[, c("rhvdr_b", "rhvdr_m", "rhvdr_t") := .(fifelse(rh_100_a0 > 5 & rh_50_a0 >=0 & rh_98_a0 >=0, rh_50_a0/rh_98_a0, -9999),
                                                         fifelse(rh_100_a0 > 5 & rh_25_a0 >=0 & rh_50_a0 >=0 & rh_75_a0 >=0 & rh_98_a0 >=0, ((rh_75_a0 - rh_25_a0)/rh_98_a0), -9999),
                                                         fifelse(rh_100_a0 > 5 & rh_50_a0 >=0 & rh_98_a0 >=0, ((rh_98_a0 - rh_50_a0)/rh_98_a0), -9999))]

  # num_modes normalized by height - num_modes_a0_hn
  # DT_m <- DT_m[, num_modes_a0_hn:=num_modes_a0/rh_98_a0]
  # PAI normalized by height - pai_a0_hn
  # DT_m <- DT_m[, pai_a0_hn:=pai_a0/rh_98_a0]
  # cover normalized by height - cover_a0_hn
  # DT_m <- DT_m[, cover_a0_hn:=cover_a0/rh_98_a0]
  # biomass normalized by height - agbd_a0_hn
  # DT_m <- DT_m[, agbd_a0_hn:=agbd_a0/rh_98_a0]

  # evenness of 1m PAI 
  # TODO: make sure RH100 > 0 (leads to a small number of Inf values)
  DT_m <- DT_m[, even_pai_1m_a0 := fhd_pai_1m_a0/log(ceiling(rh_100_a0))]

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

  DT_m[, fhd_pavd_5m_a0:= fifelse(rh_100_a0 > 5,
                                 mapply(shandiv_pavd5_func, 
                                 pavd_0_5, pavd_5_10, pavd_10_15, pavd_15_20,
                                 pavd_20_25, pavd_25_30, pavd_30_35, pavd_35_40,
                                 pavd_40_45, pavd_45_50, pavd_50_55, pavd_55_60,
                                 pavd_60_65, pavd_65_70, pavd_70_75, pavd_75_80,
                                 pavd_80_85, pavd_85_90, pavd_90_95, pavd_95_100),
                                 -9999)]
  
  DT_m[, even_pavd_5m_a0 := fifelse(rh_100_a0 > 5,
                                   mapply(evenness_pavd5_func, 
                                   fhd_pavd_5m_a0,
                                   pavd_0_5, pavd_5_10, pavd_10_15, pavd_15_20,
                                   pavd_20_25, pavd_25_30, pavd_30_35, pavd_35_40,
                                   pavd_40_45, pavd_45_50, pavd_50_55, pavd_55_60,
                                   pavd_60_65, pavd_65_70, pavd_70_75, pavd_75_80,
                                   pavd_80_85, pavd_85_90, pavd_90_95, pavd_95_100),
                                   -9999)]

  # fraction of PAVD in the bottom and top half
  DT_m <- DT_m[, pavd_bot_max_idx:=fifelse(round(rh_100_a0/2) > 5, (round((rh_100_a0/2)/5)*5)/5, 1)]
  uniq_idx = sort(unique(DT_m$pavd_bot_max_idx))
  list_DTs <- list()
  i <- 1
  for (u in uniq_idx){
    pavd_bot_cols = pavd_cols[1:u]
    DT_u <- DT_m[pavd_bot_max_idx == u]
    DT_u <- DT_u[, pavd_bot_sum:=rowSums(.SD), .SDcols = pavd_bot_cols]
    DT_u <- DT_u[, pavd_tot_sum:=rowSums(.SD), .SDcols = pavd_cols]
    DT_u <- DT_u[, pavd_bot_frac:=fifelse(rh_100_a0 > 5, pavd_bot_sum/pavd_tot_sum, -9999)]
    DT_u <- DT_u[, pavd_top_frac:=fifelse(pavd_bot_frac==-9999,-9999,1-pavd_bot_frac)]
    DT_u <- DT_u[,c("pavd_bot_sum", "pavd_tot_sum", "pavd_bot_max_idx"):=NULL]
    list_DTs[[i]] <- DT_u
    rm(DT_u)
    i <- i+1
  }
  DT_c <- rbindlist(list_DTs)
  rm(DT_m)

  # Whether or not to save a .csv
  if (!is.null(save_file)){
    # Save as csv
    fwrite(x = DT_c, row.names = FALSE, file = save_file)
    cat("Saved shots with custom metrics to", save_file, "\n")
    cat("\n")
    } else {
      # Return a data.table
      cat("Saved shots with custom metrics as data.table \n")
      cat("\n")
      return(DT_c)
    }
} # End gedi_l24a002_moremets function



# Function: gedi_shot_clip
# Author: Patrick Burns, NAU
# Last Updated: 04/22/2022
# Purpose: clip GEDI shots to a region of interest - either a polygon or box
# Arguments:
# > shot_DT = extracted L2A/B, L4A shots as either a data.table or .csv file
# > clip = a geometry to use for clipping - either a shapefile or bounding box
# > save_file = where to save the output csv. If NULL, nothing will be saved
# Returns: data.table of clipped GEDI L2A/B, L4A shots
gedi_shot_clip <- function(shot_DT = NULL,
                           clip = NULL,
                           save_file = NULL){

  # Handle shot_DT input
  type = class(shot_DT)
  if ('data.table' %in% type){ # the data.table already exists
    DT <- shot_DT
  } else if ('character' %in% type){ # read in .csv as data.table
    DT <- fread(file = shot_DT)
  }

  # Convert to combined df to spatial
  DT_sp = sp::SpatialPointsDataFrame(coords = DT[,c("lon_lm_a0", "lat_lm_a0")],
                                     data = DT,
                                     proj4string = sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
  cat("Converted GEDI shots to spatial points data frame\n")
  cat("\n")

  # Handle roi input
  # Clip to bounding box or shape
  if(tools::file_ext(clip) == "shp"){
    clip_sp = rgdal::readOGR(dsn = paste0(dirname(clip)), layer = strsplit(x = basename(clip), ".shp")[[1]])
    clip_sp_tr = sp::spTransform(x = clip_sp, CRSobj = sp::CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
    fp_clip = DT_sp[clip_sp_tr, ]
    cat("\n")
    cat("Clipped shots to this shapefile", clip, "\n")
    cat("\n")
  } else {
    inCoords = as.numeric(unlist(strsplit(clip, ",")))
    x = c(inCoords[2], inCoords[4], inCoords[4], inCoords[2], inCoords[2])
    y = c(inCoords[1], inCoords[1], inCoords[3], inCoords[3], inCoords[1])
    coords = data.frame(x=x,y=y)
    p = sp::Polygon(coords)
    ps = sp::Polygons(list(p),1)
    box = sp::SpatialPolygons(list(ps))
    sp::proj4string(box) = sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    fp_clip = DT_sp[box,]
    cat("Clipped footprints to this bounding box", inCoords, "\n")
    cat("\n")
  }

  # Check to see if there are shots after filtering
  if (nrow(fp_clip@data) == 0){
    cat("No shots within region geometry \n")
    cat("\n")
    quit(save = "no")
  } else {
    cat("Found", nrow(fp_clip@data), "shots within area of interest \n")
    cat("\n")

    # Whether or not to save a .csv
    if (!is.null(save_file)){
      # Save as csv
      fwrite(x = fp_clip@data, file = save_file)
      cat("Saved clipped shots to", save_file, "\n")

      # Save as geopackage
      # Save geopackage
      save_file_geo <- paste0(tools::file_path_sans_ext(save_file), "_geo.gpkg")
      rgdal::writeOGR(obj = fp_clip, dsn = save_file_geo, layer = "combined_sp", driver = "GPKG")
      cat("Saved clipped shots to", save_file_geo, "\n")
      cat("\n")
    } else {
      # Return a data.table
      cat("Did NOT save .csv or geopackage after clipping \n")
      cat("\n")
      return(as.data.table(fp_clip@data))
    }
  }
} # End gedi_shot_clip function



# Function: gedi_l2_join
# Author: Patrick Burns, NAU
# Last Updated: 09/28/2021
# Purpose: join gedi L2A and L2B using the shot number field
# Arguments:
# > l2a_DT = extracted L2A shots as either a data.table or .csv file
# > l2b_DT = extracted L2B shots as either a data.table or .csv file
# > join_all = whether or not to keep all rows. Default is FALSE resulting in an inner join which returns only matching shot_numbers
# > save_file = where to save the output csv. If NULL, nothing will be saved
# Returns: data.table of joined GEDI L2A+B shots
gedi_l2_join <- function(l2a_DT = NULL,
                         l2b_DT = NULL,
                         join_all = FALSE,
                         save_file = NULL){

  # Check that both l2a_DT and l2b_DT arguments are non-NULL


  # Handle l2a_DT input
  type_2a = class(l2a_DT)
  if ('data.table' %in% type_2a){ # the data.table already exists
    DT_a <- l2a_DT
  } else if ('character' %in% type_2a){ # read in .csv as data.table
    DT_a <- fread(file = l2a_DT)
  }

  # Handle l2a_DT input
  type_2b = class(l2b_DT)
  if ('data.table' %in% type_2b){ # the data.table already exists
    DT_b <- l2b_DT
  } else if ('character' %in% type_2b){ # read in .csv as data.table
    DT_b <- fread(file = l2b_DT)
  }

  # Set keys for join
  setkey(DT_a, "shot_num", "lon_lm_a0", "lat_lm_a0")
  setkey(DT_b, "shot_num", "lon_lm_a0", "lat_lm_a0")

  # Apply the join
  keep <- union(names(DT_a), names(DT_b))

  if (join_all){
    # join with all records returned
    DT_j <- DT_a[DT_b, mget(keep)]
  } else {
    # join with only matching records returned
    DT_j <- DT_a[DT_b, mget(keep), nomatch = 0]
  }

  # What to save
  # Check to see if there are shots after joining
  if (nrow(DT_j) == 0){
    cat("0 shots remain after joining \n")
    cat("\n")
    quit(save = "no")
  } else {
    cat("Found", nrow(DT_j), "L2 shots after joining L2A and L2B \n")
    cat("\n")

    # Whether or not to save a .csv
    if (!is.null(save_file)){
      # Save as csv
      fwrite(x = DT_j, row.names = FALSE, file = save_file)
      cat("Saved" , nrow(DT_j), "clipped shots to", save_file, "\n")
    } else {
      # Return a data.table
      cat("Saved", nrow(DT_j), "filtered shots as data.table \n")
      cat("\n")
      return(as.data.table(DT_j))
    }
  }
} # End gedi_l2_join function



# Function: gedi_l2l4a_join
# Author: Patrick Burns, NAU
# Last Updated: 09/26/2022
# Purpose: join gedi L2 with L4A using the shot number field
# Arguments:
# > l2_DT = extracted L2 shots as either a data.table or .csv file
# > l4a_DT = extracted L4A shots as either a data.table or .csv file
# > join_all = whether or not to keep all rows. Default is FALSE resulting in an inner join which returns only matching shot_numbers
# > save_file = where to save the output csv. If NULL, nothing will be saved
# Returns: data.table of joined GEDI L2A+B shots
gedi_l2l4a_join <- function(l2_DT = NULL,
                            l4a_DT = NULL,
                            join_all = FALSE,
                            save_file = NULL){

  # Check that both l2_DT and l4a_DT arguments are non-NULL
  # Handle l2a_DT input
  type_2 = class(l2_DT)
  if ('data.table' %in% type_2){ # the data.table already exists
    DT_2 <- l2_DT
  } else if ('character' %in% type_2){ # read in .csv as data.table
    DT_2 <- fread(file = l2_DT)
  }

  # Handle l4a_DT input
  type_4a = class(l4a_DT)
  if ('data.table' %in% type_4a){ # the data.table already exists
    DT_4a <- l4a_DT
  } else if ('character' %in% type_4a){ # read in .csv as data.table
    DT_4a <- fread(file = l4a_DT)
  }

  # Set keys for join
  setkey(DT_2, "shot_num")
  setkey(DT_4a, "shot_num")

  # Apply the join
  keep <- union(names(DT_2), names(DT_4a))

  if (join_all){
    # join with all records returned
    DT_j <- DT_2[DT_4a, mget(keep)]
  } else {
    # join with only matching records returned
    DT_j <- DT_2[DT_4a, mget(keep), nomatch = 0]
  }

  # What to save
  # Check to see if there are shots after joining
  if (nrow(DT_j) == 0){
    cat("0 shots remain after joining \n")
    cat("\n")
    quit(save = "no")
  } else {
    cat("Found", nrow(DT_j), "shots after joining L2 and L4A \n")
    cat("\n")

    # Whether or not to save a .csv
    if (!is.null(save_file)){
      # Save as csv
      fwrite(x = DT_j, row.names = FALSE, file = save_file)
      cat("Saved" , nrow(DT_j), "clipped shots to", save_file, "\n")
    } else {
      # Return a data.table
      cat("Saved", nrow(DT_j), "filtered shots as data.table \n")
      cat("\n")
      return(as.data.table(DT_j))
    }
  }
} # End gedi_l2l4a_join function



# Function: gedi_outliers_umd
# Author: Patrick Burns, NAU
# Last Updated: 10/08/2022
# Purpose: remove outlier granules based on UMD local outlier detection procedure
# Arguments:
# > table = data.table or txt/csv file with GEDI shots
# > grid_shp = the shapefile grid that is associated with outliers
# > grid_name = the name of the field used to identify each grid
# > outliers_tab = a .csv file that lists bad orbit-grid combinations
# > remove = whether or not to actually remove the outliers (TRUE/FALSE)
# > save_file = where to save the output csv. If NULL, nothing will be saved
# Returns: data.table of quality shots
gedi_outliers_umd <- function(table = NULL,
                               grid_shp = '/projects/above_gedi/users/pburns/GEDI/grids/grid_EASE2_72km_land_pm52_named.shp',
                               grid_name = 'EASE72_id',
                               outliers_tab = '/projects/above_gedi/users/pburns/GEDI/bad_granules/issgedi_l4b_excluded_granules_r002_thruMW222_rec20230917.csv',
                               remove = FALSE,
                               save_file = NULL){
  
  # Handle DT input
  type = class(table)
  if ('data.table' %in% type){ # the data.table already exists
    DT <- table
  } else if ('character' %in% type){ # read in .csv as data.table
    DT <- fread(file = table)}
  
  # Load the grid shapefile
  grid <- st_read(grid_shp)
  
  # Load the outlier table
  DT_outliers <- fread(outliers_tab)
  outliers_vec <- DT_outliers[[1]]
  
  # GEDI shots as sf and transformed to CRS of grid
  DT_sf <- st_as_sf(DT, coords = c("lon_lm_a0", "lat_lm_a0"),
                    crs = 'EPSG:4326', remove = FALSE) %>%
    st_transform(st_crs(grid))
  
  # Spatially (inner) join the GEDI shots with grid
  DT_sf_j <- as.data.table(st_join(DT_sf, grid, left = FALSE))
  
  # Split the geometry field into lon and lat, then remove geometry
  DT_sf_j <- DT_sf_j[ , lon_lm_a0_6933 := sf::st_coordinates(geometry)[,1]]
  DT_sf_j <- DT_sf_j[ , lat_lm_a0_6933 := sf::st_coordinates(geometry)[,2]]
  DT_sf_j <- DT_sf_j[ , geometry := NULL]
  
  # add a field for excluding orbit-grid outliers and then exclude them
  DT_sf_j <- DT_sf_j[, orbit_num := as.numeric(substr(orbit,2,6))]
  DT_sf_j <- DT_sf_j[, ease72_og := paste0(orbit_num, sub_orbit, "-", get(grid_name))]
  DT_sf_j <- DT_sf_j[, loc_out_umd  := ifelse(ease72_og %in% outliers_vec,1,0)]
  DT_sf_j <- DT_sf_j[, EASE72_id := NULL]
  DT_sf_j <- DT_sf_j[, orbit := NULL]
  DT_sf_j <- DT_sf_j[, orbit := orbit_num]
  DT_sf_j <- DT_sf_j[, orbit_num := NULL]
  nrow_orig <- nrow(DT_sf_j)
  if (remove){
    DT_sf_j_o <- DT_sf_j[loc_out_umd == 0]
    nrow_noout <- nrow(DT_sf_j_o)
  } else{
    DT_sf_j_o <- DT_sf_j
  }
  cat("Identified", nrow(DT_sf_j[loc_out_umd == 1]), "outlier shots \n")
  
  
  # What to save
  # Check to see if there are shots after joining
  if (nrow(DT_sf_j_o) == 0){
    cat("0 shots remain after excluding bad granules \n")
    cat("\n")
    quit(save = "no")
  } else {
    if (remove){
      cat("Removed outliers -",nrow_noout, "shots remain \n")  
    } else {
      cat("Did not remove outliers -",nrow_orig, "shots remain \n")
    }
    cat("\n")

    # Whether or not to save a .csv
    if (!is.null(save_file)){
      # Save as csv
      fwrite(x = DT_sf_j_o, row.names = FALSE, file = save_file)
      cat("Saved shots to", save_file, "\n")
    } else {
      # Return a data.table
      cat("Saved shots as data.table \n")
      cat("\n")
      return(as.data.table(DT_sf_j_o))
    }
  }
} # End gedi_outliers_umd function



# Function: gedi_l2l4a_hqf
# Author: Patrick Burns, NAU
# Last Updated: 10/14/2023
# Purpose: add a high quality flag for vegetation and optionally remove low quality shots
# Arguments:
# > table = data.table or txt/csv file with GEDI shots
# > remove = whether or not to actually remove the poor quality shots (TRUE/FALSE)
# > save_file = where to save the output csv. If NULL, nothing will be saved
# Returns: a quality-filtered data.table or 
gedi_l2l4a_hqf <- function(table = NULL,
                           remove = FALSE,
                           save_file = NULL){
  # Handle DT input
  type = class(table)
  if ('data.table' %in% type){ # the data.table already exists
    DT <- table
  } else if ('character' %in% type){ # read in .csv as data.table
    DT <- fread(file = table)}
  
  # Apply these filters
  # > degrade flag (accept these values)
  deg_ok <- c(0,3,8,10,13,18,20,23,28,30,33,38,40,43,48,60,63,68)
  # > not leaf off
  # > surface water percentage less than 10%
  # > urban proportion less than 50%
  # valid values of pai and cover
  DT_hqf <- DT[, l2_hqflag:=ifelse(deg_flag %in% deg_ok 
                                    & l2b_algrun_flag == 1
                                    & l2b_qflag_a0 == 1
                                    & ls_waterp < 10
                                    & urb_prop < 50 
                                    & !(leafoff_flag == 1)
                                    & loc_out_umd == 0
                                    & pai_a0 >= 0
                                    & pai_l1 >= 0
                                    & pavd_0_5 >= 0
                                    & cover_a0 >= 0
                                    & cover_a0 <= 1
                                    & cover_l1 >= 0
                                    & cover_l1 <= 1
                                    & abs(elev_diff_dem) < 150, 
                                    1, 0)]

  DT_hqf <- DT[, l4a_hqflag:=ifelse(l2_hqflag == 1
                                    & l4a_algrun_flag == 1
                                    & l4a_qflag_a0 == 1
                                    & agbd_a0 >= 0,
                                    1, 0)]
  nrow_orig <- nrow(DT_hqf)

  if (remove){
    DT_o <- DT_hqf[l2_hqflag == 1]
    nrow_noout <- nrow(DT_o)
  } else{
    DT_o <- DT_hqf
  }

  # What to save
  # Check to see if there are shots after joining
  if (nrow(DT_o) == 0){
    cat("0 shots remain after excluding bad granules \n")
    cat("\n")
    quit(save = "no")
  } else {
    if (remove){
      cat("Removed low quality shots -",nrow_noout, "shots remain \n")  
    } else {
      cat("Did not remove low quality shots -",nrow_orig, "shots remain \n")
    }
    cat("\n")

    # Whether or not to save a .csv
    if (!is.null(save_file)){
      # Save as csv
      fwrite(x = DT_o, row.names = FALSE, file = save_file)
      cat("Saved shots to", save_file, "\n")
    } else {
      # Return a data.table
      cat("Saved shots as data.table \n")
      cat("\n")
      return(as.data.table(DT_o))
    }
  }
} # End gedi_l2l4a_hqf function



# Function: gedi_chunk_gcs
# Author: Patrick Burns, NAU
# Last Updated: 10/18/2022
# Purpose: chunk GEDI shots using a regularly spaced grid
# Arguments:
# > table = GEDI shots as either a data.table or .csv file
# > chunk_dim = square dimension in decimal degrees to use for splitting
# > chunk_x_max = Maximum longitude. Will be mirrored to -x_max
# > chunk_y_max = maximum latitude. Will be mirrored to -y_max
# > save_file = where to save the output csv. If NULL, nothing will be saved
# > summ_file = where to save the summary file
# Returns: data.table of clipped GEDI L2 (A and/or B) shots

gedi_chunk_gcs <- function(table = NULL,
                           chunk_dim = 1.0,
                           chunk_x_max = 180.0,
                           chunk_y_max = 52.0,
                           save_file = NULL,
                           summ_file = NULL){

  # Handle l2_DT input
  type = class(table)
  if ('data.table' %in% type){ # the data.table already exists
    DT <- table
  } else if ('character' %in% type){ # read in .csv as data.table
    DT <- fread(file = table)
  }

  if (nrow(DT) == 0){
    cat("0 shots for gridding \n")
    cat("\n")
    quit(save = "no")
  }

  # Write csvs for specified grid
  xseq <- seq(chunk_x_max*-1, (chunk_x_max-chunk_dim), chunk_dim)
  yseq <- seq(chunk_y_max*-1, (chunk_y_max-chunk_dim), chunk_dim)

  # Limit the x and y sequences by min and max GEDI coordinates
  gedi_x_min <- floor(min(DT$lon_lm_a0, na.rm = TRUE))
  gedi_x_max <- floor(max(DT$lon_lm_a0, na.rm = TRUE))
  gedi_y_min <- floor(min(DT$lat_lm_a0, na.rm = TRUE))
  gedi_y_max <- floor(max(DT$lat_lm_a0, na.rm = TRUE))

  xseq_f <- xseq[xseq >= gedi_x_min & xseq <= gedi_x_max]
  yseq_f <- yseq[yseq >= gedi_y_min & yseq <= gedi_y_max]


  # Loop through XY Grids
  for (xmin in xseq_f){
    xmax <- as.numeric(xmin + chunk_dim)
    if (xmin < 0){
      EW <- "W"
    } else {
      EW <- "E"
    }

    for (ymin in yseq_f){
      ymax <- as.numeric(ymin + chunk_dim)
      if (ymin < 0){
        NS <- "S"
      } else {
        NS <- "N"
      }

      # Construct the tile ID (corresponds to lower left corner)
      EW_pad <- paste0(sprintf("%03d", abs(xmin)), EW)
      NS_pad <- paste0(sprintf("%02d", abs(ymin)), NS)
      tileid <- paste0(EW_pad, NS_pad)

      # Filter GEDI shots by bounds
      DT_g <- DT[lon_lm_a0 >= xmin & lon_lm_a0 < xmax & lat_lm_a0 >= ymin & lat_lm_a0 < ymax]

      if (nrow(DT_g) > 0){
        cat(tileid, "-", nrow(DT_g), "shots \n")
        # Save as csv
        chunk_save_file <- paste0(tools::file_path_sans_ext(save_file), "_g", tileid, ".csv")
        fwrite(x = DT_g, row.names = FALSE, file = chunk_save_file)

        if (!is.null(summ_file)){
          DT_s <- data.table(file = paste0(tools::file_path_sans_ext(save_file), "_g", tileid, ".csv"),
                             n_filt_shots = nrow(DT_g),
                             n_l2_hq_shots = nrow(DT_g[l2_hqflag == 1]),
                             n_l4a_hq_shots = nrow(DT_g[l4a_hqflag == 1]),
                             n_locout_shots = nrow(DT_g[loc_out_umd == 1]),
                             xmin_shots = min(DT_g$lon_lm_a0, na.rm = TRUE),
                             xmax_shots = max(DT_g$lon_lm_a0, na.rm = TRUE),
                             ymin_shots = min(DT_g$lat_lm_a0, na.rm = TRUE),
                             ymax_shots = max(DT_g$lat_lm_a0, na.rm = TRUE),
                             date_dec_min = min(DT_g$date_dec, na.rm = TRUE),
                             date_dec_max = max(DT_g$date_dec, na.rm = TRUE)
                             )
          chunk_summ_save_file <- paste0(tools::file_path_sans_ext(summ_file), "_g", tileid, ".csv")
          fwrite(x = DT_s, row.names = FALSE, file = chunk_summ_save_file)
        }
        #cat("\n")
        
      }
      rm(DT_g)
      
    } # end Ys
    
  } # end Xs and chunking loop
  cat("\n")
  
} # End gedi_chunk_gcs function



# Function: gedi_comb_chunk_summs
# Author: Patrick Burns, NAU
# Last Updated: 10/10/2022
# Purpose: write file list of processed GEDI shots that are associated with unique chunks
# Arguments:
# > summ_comb_tab = the combined .csv file
# > chunk_add_dist = the distance from the main chunk to search for neighboring files
# > save_dir = where to save output neighboring file lists
# Returns: .csv file(s) of neighboring files
gedi_comb_chunk_summs <- function(summ_comb_tab = NULL,
                                  chunk_add_dist = 0.5,
                                  save_dir = NULL){

  # read in the combined summary table
  DT_s <- fread(summ_comb_tab)

  # Get the unique chunk names
  chunks <- sort(unique(tools::file_path_sans_ext(x = tstrsplit(x = DT_s[[1]], split = "_filt_", fixed = TRUE)[[2]])))

  # loop over each chunk to find neighbors with a certain distance of the chunk of interest
  for (c in chunks){
    cat("Working on chunk", c, "\n")
    
    xs = substr(c, 5, 5)
    if (xs == "W"){
      xmin <- (-1*as.numeric(substr(c, 2, 4)))
    } else if(xs == "E"){
      xmin <- as.numeric(substr(c, 2, 4))
    }
    xmin_wind <- as.numeric(xmin-chunk_add_dist)
    xmax_wind <- as.numeric(xmin+1+chunk_add_dist)

    ys = substr(c, 8, 8)
    if (ys == "S"){
      ymin <- (-1*as.numeric(substr(c, 6, 7)))
    } else if (ys == "N"){
      ymin <- as.numeric(substr(c, 6, 7))
    }
    ymin_wind <- as.numeric(ymin-chunk_add_dist)
    ymax_wind <- as.numeric(ymin+1+chunk_add_dist)
    
    # Select the files that overlap the specified window and save the GEDI shot file names associated with each chunk
    # OLD (doesnt work): DT_s[xmin_shots >= xmin_wind & xmax_shots <= xmax_wind & ymin_shots >= ymin_wind & ymax_shots <= ymax_wind]
    DT_s_w <- DT_s[xmax_shots >= xmin_wind & xmin_shots <= xmax_wind & ymax_shots >= ymin_wind & ymin_shots <= ymax_wind]
    fwrite(DT_s_w[,1], paste0(save_dir, c, "_ngb_files.csv"), col.names = FALSE)
  }
} # End gedi_comb_chunk_summs function