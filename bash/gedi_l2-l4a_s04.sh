#!/bin/bash

# ----- ABOUT -----
# Title: gedil2l4_s04.sh
# Author: Patrick Burns [pb463@nau.edu]
# Purpose: mosaic GEDI rasterized chunks
# Last Updated: 24 Oct. 2022
# TODO:



# ----- INPUTS -----
# the main processing directory (should contain s03)
procDir='/projects/geode_data/share/lidar/GEDI0204a.002/GEDIv002_L0204A_20190417to20230316_proc202310/'

# date pairs to process. First and second date should be separate by _. Pairs should be separated by commas.
date_pairs='20190417_20191231,20200101_20201231,20210101_20211231,20220101_20221231,20230101_20230316'

# CRS EPSG (4326 or 6933)
crs_epsg=6933


## SLURM settings
# Processing time to request for script 1
# try 6 to 8 hrs for global
timeReq_s04="08:00:00"

# Memory to request for script 3
# 25G should be enough, might need up to 40G
memoryReq_s04="40G"

# Specify how many array jobs to try and do at once (200 seems to work ok)
#
arrayLim=50



# ----- PROCESSING -----
echo

# Starting date and time of the script for tagging log names
date_time=`date +%Y%m%d_%H%M%S`

# the directory tree for s04
s04Dir="${procDir}s04/"
if [ ! -d "$s04Dir" ]
then
    mkdir "$s04Dir"
    mkdir "${s04Dir}vrt/"
    mkdir "${s04Dir}tif/"
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


# make the list of metric-date combinations to process
fList="${s04Dir}metrics_dates_flist.csv"

metrics=$(find "${procDir}s03/" -maxdepth 2 -type f -name 'gediv002'*'.tif' -printf '%f\n' | cut -d "_" -f2,3 | sort -T /scratch/pb463/ | uniq)
#date1=$(find $procDir's03/' -maxdepth 2 -type f -name *'stats'*$crs_epsg'.tif' -printf '%f\n' | cut -d "_" -f4 | sort | uniq)
#date2=$(find $procDir's03/' -maxdepth 2 -type f -name *'stats'*$crs_epsg'.tif' -printf '%f\n' | cut -d "_" -f5 | sort | uniq)
#date_pairs=$(for d1 in $date1; do for d2 in $date2; do echo $d1'_'$d2; done; done)
IFS=', ' read -r -a date_pairs_arr <<< "$date_pairs"
res=$(find $procDir's03/' -maxdepth 2 -type f -name 'gediv002'*'.tif' -printf '%f\n' | cut -d "_" -f6 | cut -d "m" -f1 | sort -T /scratch/pb463/ | uniq | sort -T /scratch/pb463/ -r | tr "\n" " ")

# remove old file list
if [ -f "$fList" ]
then
  rm "$fList"
fi

# make new file list
touch "$fList"
for m in $metrics
do
  for d in "${date_pairs_arr[@]}"
  do
    echo "${m},${d}" 
  done
done > "$fList"

# Count files for SLURM Array
numFiles=$(wc -l "$fList" | cut -d " " -f1)


# ----- SLURM ARRAY -----
# Create the slurm array job script in the slurm folder
cd "${slurmDir}scripts/"
script_s04="sbatch_s04_${date_time}.sh"

# Specify the basename for the output log file
slurmLog_s04="${slurmDir}logs/s04_${date_time}_%A_%a.out"

cat > "$script_s04" <<EOT
#!/bin/bash

# - SLURM -
#SBATCH --job-name=gedil2l4_s04
#SBATCH --output="$slurmLog_s04"
#SBATCH --time="$timeReq_s04"
#SBATCH --mem="$memoryReq_s04"
#SBATCH --array=1-$((numFiles))
#'%'$arrayLim

SECONDS=0
date_time=`date +%Y%m%d_%H%M%S`
echo "The starting date_time: " \$date_time
echo
#echo "SLURM_JOBID: "\$SLURM_JOBID
echo "SLURM_ARRAY_JOB_ID: "\$SLURM_ARRAY_JOB_ID
echo "SLURM ARRAY TASK ID: "\$SLURM_ARRAY_TASK_ID
echo

# Activate the necessary modules and environments on Monsoon
module load anaconda3
conda activate /projects/above_gedi/users/pburns/envs/geospatial

# Get the metric and date pair for making the mosaic
line=\$(sed \$SLURM_ARRAY_TASK_ID'q;d' "$fList")
met=\$(echo \$line | cut -d "," -f1)
date_pair=\$(echo \$line | cut -d "," -f2)

echo "Looping over these spatial resolutions: ${res}"
echo
for r in $res
do
  echo "Mosaicing "\$met" for date range "\$date_pair" at a resolution of "\$r"m"
  mos_base='gediv002_'\$met'_'\$date_pair'_'\$r'm'
  # build the vrt
  cd $procDir
  srun gdalbuildvrt $s04Dir'vrt/'\$mos_base'.vrt' 's03/'*'/'\$mos_base'.tif' -overwrite -srcnodata 'nan' -vrtnodata -9999

  # rename bands in the VRT
  i=1
  if [ "\$met" = "counts_ga" ] || [ "\$met" = "counts_va" ]
  then
    for b in shots_count orbits_uniq tracks_uniq shots_nni
    do
      search="band=\"\${i}\""
      replace="band=\"\${i}\" Description=\"\${b}\""
      sed -i -e "s/\$search/\$replace/g" $s04Dir'vrt/'\$mos_base'.vrt'
      i=\$(( i + 1 ))
    done
  else
    for b in mean meanbse median sd iqr p95 shan countf
    do 
      search="band=\"\${i}\""
      replace="band=\"\${i}\" Description=\"\${b}\""
      sed -i -e "s/\$search/\$replace/g" $s04Dir'vrt/'\$mos_base'.vrt'
      i=\$(( i + 1 ))
    done
  fi

  # save the final tif (remove the older version if necessary)
  if [ -f $s04Dir'tif/'\$mos_base'.tif' ]
  then
    rm $s04Dir'tif/'\$mos_base'.tif'
  fi
  srun gdal_translate -stats -of cog -co BLOCKSIZE=256 -co COMPRESS=LZW -co BIGTIFF=IF_SAFER -co OVERVIEW_RESAMPLING=NEAREST $s04Dir'vrt/'\$mos_base'.vrt' $s04Dir'tif/'\$mos_base'.tif'
  
  # replace nans
  #srun gdal_calc.py -A $s04Dir'tif/'\$mos_base'_0.tif' --outfile=$s04Dir'tif/'\$mos_base'.tif' --calc="nan_to_num(A, nan=-9999)" --NoDataValue=-9999
  # rm $s04Dir'tif/'\$mos_base'_0.tif'
  echo
done


# - ENDING -
date_time=`date +%Y%m%d_%H%M%S`
echo "The ending date_time: " \$date_time
duration=\$SECONDS
echo "\$((\$duration / 60)) minutes and \$((\$duration % 60)) seconds elapsed."

EOT

# Run the slurm array job script
echo "----- Starting Script 4: Mosaic chunks into one raster -----"
jobid_s04=$(sbatch --parsable "$script_s04")
echo "The SLURM job id is: ${jobid_s04}"
