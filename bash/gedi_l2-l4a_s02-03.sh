#!/bin/bash


# ----- ABOUT -----
# Title: gedil2l4a_s02-03.sh
# Author: Patrick Burns [pb463@nau.edu]
# Purpose: combine GEDI processed quality shot tables using a regular grid and rasterize select metrics
# Last Updated: 10 July 2023
# TODO:



# ----- INPUTS -----
# The base processing directory
procDir='/scratch/pb463/gedi/L2A-L2B-L4A_002_2019-04-17_2023-03-16_global/'

# an alternative path to save the results (optional).
# Leave blank if results of s03 should be saved in the processing directory
alt_saveDir='/projects/geode_data/share/lidar/GEDI0204a.002/GEDIv002_L0204A_20190417to20230316_proc202310/'

# GEDI metrics to process (comma separated)
gedi_metrics='date_dec,sens_a0,elev_lm_a0,num_modes_a0,rh_50_a0,rh_95_a0,rh_98_a0,pai_a0,fhd_pai_1m_a0,cover_a0,agbd_a0,pavd_0_5_frac,pavd_max_h,rhvdr_b,rhvdr_m,rhvdr_t,pavd_bot_frac,pavd_top_frac,fhd_pavd_5m_a0,even_pai_1m_a0,even_pavd_5m_a0,pavd_0_5,pavd_5_10,pavd_10_15,pavd_15_20,pavd_20_25,pavd_25_30,pavd_30_35,pavd_35_40,pavd_40_45,pavd_45_50,pavd_50_55,pavd_55_60,pavd_60_65,pavd_65_70,pavd_70_75,pavd_75_80'

# the global base raster
base_raster='/scratch/pb463/EASE/EASE2_000300m_01_int16.tif'
#'/scratch/pb463/FLII/flii_earth_byte.tif'

# array of additional spatial resolutions to resample the global base raster to
# should be in the units of the base raster (either decimal degrees or meters)
res_arr=( 1000 6000 12000 )
#( 0.008983112 0.05389867 0.107797341 )

# distance (in decimal degrees) to search from each chunk (should be about 2x the maximum pixel size)
chunk_add_dist=0.25

# date pairs
date_pairs='2019-04-17,2019-12-31,2020-01-01,2020-12-31,2021-01-01,2021-12-31,2022-01-01,2022-12-31,2023-01-01,2023-03-16'
#'2019-04-17,2023-03-16'
#

## SLURM settings
# Processing time to request for script 1
# 07:00:00 works well for 30 deg. lat and less
# 06:00:00 works well for >30 deg. lat (full mission)
timeReq_s03='07:00:00'

# Memory to request for script 3
# 50G works well for 30 deg. lat and less
# 60G works well for annual
# 80G works well for >30 deg. lat (full mission)
memoryReq_s03='60G'

# Specify how many array jobs to try and do at once (200 seems to work ok)
arrayLim=500



# ----- PROCESSING -----
echo

# Starting date and time of the script for tagging log names
date_time=`date +%Y%m%d_%H%M%S`

# the path containing the orbital summaries
summ_dir="${procDir}s01/summaries/"

# the directory tree for s02
s02Dir="${procDir}s02/"
if [ ! -d "$s02Dir" ]
then
    mkdir "$s02Dir"
    mkdir "${s02Dir}ngb_lists/"
fi

# whether or not to save the s03 results in the current processing directory or somewhere else
if [ ! -z "$alt_saveDir" ]
then
  s03Dir="${alt_saveDir}s03/"
  if [ ! -d "$alt_saveDir" ]
  then
    mkdir -p "$alt_saveDir"
  fi
else
  s03Dir="${procDir}s03/"
fi

if [ ! -d "$s03Dir" ]
then
    mkdir "$s03Dir"
    mkdir "${s03Dir}base_rasters/"
fi

# Make slurm subdirectories
slurmDir="${procDir}slurm/"
if [ ! -d "$slurmDir" ]
then
    mkdir "$slurmDir"
    mkdir "${slurmDir}scripts/"
    mkdir "${slurmDir}jobstats/"
    mkdir "${slurmDir}logs/"
fi

# slurm subdirectories
slurmDir="${procDir}slurm/"
logDir="${slurmDir}logs/"
statsDir="${slurmDir}jobstats/"

module load anaconda3
conda activate /projects/above_gedi/users/pburns/envs/geospatial

# Run R script to identifying neighboring chunks
fList="${s02Dir}ngb_flist_full.csv"
if [ ! -f "$fList" ]
then
  echo "----- Starting Script 2: Find neighbors -----"
  summ_tab_comb="${s02Dir}summ_tabs_comb.csv"
  if [ ! -f "$summ_tab_comb" ]
  then
    echo "Merging all summary .csv files into one table (can be slow for millions of files): ${summ_tab_comb}"
    # add the header to the merged file
    echo 'file,n_filt_shots,n_l2_hq_shots,n_l4a_hq_shots,n_locout_shots,xmin_shots,xmax_shots,ymin_shots,ymax_shots,date_dec_min,date_dec_max' > "$summ_tab_comb"
    # find all summary files and take the 2nd line
    find "$summ_dir" -type f -name *.csv -exec tail -q -n +2 {} \; >> "$summ_tab_comb"
    echo "Finished merging all summary .csv files into one table: ${summ_tab_comb}"
    echo
  else
    echo "Merged summary table already exists: ${summ_tab_comb}"
    echo
  fi

  slurmLog_s02a="${slurmDir}logs/s02a_${date_time}_%A.out"
  echo "Starting script to find neighboring tables..."
  srun --mem=10G --time=02:00:00 --output="$slurmLog_s02a" Rscript /home/pb463/scripts/repos/gedi_gridding/R/gedi_l2-l4a_s02_find_ngbs.R "$summ_tab_comb" "$chunk_add_dist" "${s02Dir}ngb_lists/"
  echo "Finished script to find neighboring tables."
  echo

  # make a list of the files that list neighbors for each grid
  find "${s02Dir}ngb_lists/" -type f -name "*.csv" > "$fList"
  sort -o "$fList" "$fList"
fi

# Count files for SLURM Array
numFiles=$(wc -l "$fList" | cut -d " " -f1)

# prep global base rasters
echo "----- Prepping global base rasters for Script 3 -----"
# copy the original raster
cp "$base_raster" "${s03Dir}base_rasters/"$(basename "$base_raster")
echo "Copied original global base raster"
for res in "${res_arr[@]}"
do
  file_base=$(basename "$base_raster" | cut -d "." -f1)
  base_raster_resamp="${s03Dir}base_rasters/${file_base}_${res}.tif"
  gdal_translate -ot Byte -tr "$res" "$res" -r nearest -q -a_nodata none "$base_raster" "$base_raster_resamp" -co COMPRESS=LZW -co TILED=YES
  echo "Created another global base raster at a spatial resolution of ${res}"
done
echo
# remove the finest resolution raster
rm "${s03Dir}base_rasters/"$(basename "$base_raster")


# ----- SLURM ARRAY -----
# Create the slurm array job script in the slurm folder
cd "${slurmDir}scripts/"
script_s03="sbatch_s03_${date_time}.sh"

# Specify the basename for the output log file
slurmLog_s03="${slurmDir}logs/s03_${date_time}_%A_%a.out"

cat > "$script_s03" <<EOT
#!/bin/bash

# - SLURM -
#SBATCH --job-name=gedil2l4_s03
#SBATCH --output="$slurmLog_s03"
#SBATCH --time="$timeReq_s03"
#SBATCH --mem="$memoryReq_s03"
#SBATCH --array=1-$((numFiles))'%'"$arrayLim"


SECONDS=0
date_time=`date +%Y%m%d_%H%M%S`
echo "The starting date_time: "\$date_time
echo
#echo "SLURM_JOBID: "\$SLURM_JOBID
echo "SLURM_ARRAY_JOB_ID: "\$SLURM_ARRAY_JOB_ID
echo "SLURM ARRAY TASK ID: "\$SLURM_ARRAY_TASK_ID
echo

# Activate the necessary modules and environments on Monsoon
module load anaconda3
conda activate /projects/above_gedi/users/pburns/envs/geospatial

# Get the file that lists all touching GEDI orbits for a particular chunk
inFile=\$(sed \$SLURM_ARRAY_TASK_ID'q;d' $fList)

# Make a sub-folder to process in
chunk_id=\$(basename \$inFile | cut -d "_" -f1)
chunk_id_low=\$(echo \$chunk_id | tr '[:upper:]' '[:lower:]')
echo "Chunk ID: "\$chunk_id
chunk_dir=$s03Dir\$chunk_id'/'
if [ ! -d \$chunk_dir ]
then
    mkdir \$chunk_dir
fi

# re-make the base_rasters dir in the event that jobs fail
if [ -d \$chunk_dir'base_rasters/' ]
then
  rm -r \$chunk_dir'base_rasters/'
  mkdir \$chunk_dir'base_rasters/'
else
  mkdir \$chunk_dir'base_rasters/'
fi

# Get coords for the bounding box and buffer
EW=\$(echo \$chunk_id | cut -c5)
if [ \$EW = "W" ]
then
    EW_mult=-1
else
    EW_mult=1
fi

NS=\$(echo \$chunk_id | cut -c8)
if [ \$NS = "S" ]
then
    NS_mult=-1
else
    NS_mult=1
fi

xmin_str_ld=\$(echo \$chunk_id | cut -c4)
xmin_str=\$(echo \$chunk_id | cut -c2-3 | sed 's/^0*//')\$xmin_str_ld
xmin=\$(( \$xmin_str * \$EW_mult ))
xmin_w=\$(echo \$xmin - $chunk_add_dist | bc)
xmax_w=\$(echo \$xmin + 1.0 + $chunk_add_dist | bc)

ymin_str_ld=\$(echo \$chunk_id | cut -c7)
ymin_str=\$(echo \$chunk_id | cut -c6 | sed 's/^0*//')\$ymin_str_ld
ymin=\$(( \$ymin_str * \$NS_mult ))
ymin_w=\$(echo \$ymin - $chunk_add_dist | bc)
ymax_w=\$(echo \$ymin + 1.0 + $chunk_add_dist | bc)
echo "Buffered bounding box: "\$xmin_w \$xmax_w \$ymin_w \$ymax_w
echo


# combine .csvs for the particular chunk
tab_comb=\$chunk_dir'tab_ngb_comb_'\$chunk_id'.csv'
tab1=\$(sed '1q;d' \$inFile)
head -n 1 \$tab1 > \$tab_comb
tail -n +2 -q \$(cat \$inFile) >> \$tab_comb
n_tabs=\$(wc -l \$inFile | cut -d " " -f1)
n_lines=\$(wc -l \$tab_comb | cut -d " " -f1)
echo "Merged "\$n_tabs" neighboring .csv files into one table:"
echo \$tab_comb
echo "with "\$n_lines" total shots"
echo


# crop the global base rasters to the buffered chunk
echo "Preparing rasters for rasterization..."
bases=($s03Dir'base_rasters/'*.tif)
for ras in \${bases[@]}
do
  file_base=\$(basename \$ras)
  base_crop=\$chunk_dir'base_rasters/'\$file_base
  gdal_translate -a_nodata none -projwin \$xmin_w \$ymax_w \$xmax_w \$ymin_w -projwin_srs 'EPSG:4326' \$ras \$base_crop -co COMPRESS=LZW -co TILED=YES
  echo "Cropped base raster "\$ras" to use for rasterization"
  echo
done
echo

# rasterize the GEDI shots and limit the final table to the chunk bounds
echo "Starting R rasterization script..."
Rscript ~/scripts/repos/gedi_gridding/R/gedi_l2-l4a_s03_rasterize_chunks.R \$chunk_id \$tab_comb "$gedi_metrics" "$date_pairs" \$chunk_dir'base_rasters/' "$chunk_add_dist" \$chunk_dir
echo

# zip the table corresponding to the 1x1 deg. bounds
tab_comb_r=\$chunk_dir'gediv002_l2l4a_va_'\$chunk_id_low'.csv'
if [ -f \$tab_comb_r ]
then
  tab_comb_r_zip=\$chunk_dir'gediv002_l2l4a_va_'\$chunk_id_low'.zip'
  zip -q -j \$tab_comb_r_zip \$tab_comb_r
  rm \$tab_comb_r
fi
echo 'Zipped '\$tab_comb_r
#rm \$tab_comb_r_zip
echo

rm -r \$chunk_dir'base_rasters/'
rm \$tab_comb
echo 'Removed extra files'
echo

# - ENDING -
date_time=`date +%Y%m%d_%H%M%S`
echo "The ending date_time: "\$date_time
duration=\$SECONDS
echo "\$((\$duration / 60)) minutes and \$((\$duration % 60)) seconds elapsed."

EOT

# Run the slurm array job script
echo "----- Starting Script 3: Merge neighboring tables and rasterize -----"
jobid_s03=$(sbatch --parsable "$script_s03")
echo "The SLURM job id is: ${jobid_s03}"
