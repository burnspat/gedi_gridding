#!/bin/bash

# ----- ABOUT -----
# Title: gedil2l4_s00-01.sh
# Author: Patrick Burns [pb463@nau.edu]
# Purpose: Uses SLURM arrays to extract, filter, join, and grid GEDI L2 and L4a data
# Last Updated: 14 Oct. 2023
# TODO:



# ----- INPUTS -----
# GEDI products to process. Options are:
# L2A, L2B, L4A, L2A-L2B, L2A-L2B-L4A
prod="L2A-L2B-L4A"

# directories that GEDI data were downloaded to
l2a_dir="/projects/geode_data/share/lidar/GEDI02_A.002/0_orig/alldata/"
l2b_dir="/projects/geode_data/share/lidar/GEDI02_B.002/0_orig/alldata/"
l4a_dir="/projects/geode_data/share/lidar/GEDI04_A.002_1/0_orig/data/"

# Unique name for the region (no spaces or underscores)
region_name='global'

# Shapefile or bounding box to use for clipping GEDI shots
region_geom="52,-180,-52,180"

# The spatial dimension to use for gridding. Use 'NULL' to skip gridding
chunk_dim=1.0

# Start date (YYYY-MM-DD)
s_date='2019-04-17'

# End date (YYYY-MM-DD)
e_date='2023-03-16'

# Minimum day of year
min_doy=1

# Maximum day of year
max_doy=366

# Processing directory (this will be created by the script if it doesn't already exist)
baseDir='/scratch/pb463/gedi/'

# Vector for density summary (leave NULL if not needed)
dens_vec='NULL'

# Vector attribute for density summary (leave NULL if not needed)
vec_attr='NULL'


## SLURM settings
# Processing time to request for script 1
# usually takes less than 10 min for a few tasks, but need to use at least 20 min when the array limit is higher (>200)
timeReq_s1="00:20:00"

# Memory to request for script 1
# 12G is usually a good starting point for L2AB + L4A
memoryReq_s1="12G"

# Specify how many array jobs to try and do at once (200 seems to work ok)
# 300 seems ok (1% failure rate). 500 works, but a higher percentage of jobs fail (~10%)
arrayLim=300



# ----- PROCESSING -----

# Starting date and time of the script for tagging log names
date_time=`date +%Y%m%d_%H%M%S`

## Make the folder structure
if [ ! -d "$baseDir" ]
then
  mkdir "$baseDir"
fi

procID="${prod}_002_${s_date}_${e_date}_${region_name}"
procDir="${baseDir}${procID}"/
if [ ! -d "$procDir" ]
then
    mkdir "$procDir"
fi

# Make subdirectories
s00Dir="${procDir}s00"/
if [ ! -d "$s00Dir" ]
then
    mkdir "$s00Dir"
fi

s01Dir="${procDir}s01"/
if [ ! -d "$s01Dir" ]
then
    mkdir "$s01Dir"
    mkdir "${s01Dir}tables"/
    mkdir "${s01Dir}summaries"/
fi

# Make slurm subdirectories
slurmDir="${procDir}slurm"/
if [ ! -d "$slurmDir" ]
then
    mkdir "$slurmDir"
    mkdir "${slurmDir}scripts"/
    mkdir "${slurmDir}jobstats"/
    mkdir "${slurmDir}logs"/
fi

# Specify the output log folder
logDir="${slurmDir}logs"/
statsDir="${slurmDir}jobstats"/


# Run R script to spatially query GEDI granules and filter by date
fList="${s00Dir}${procID}_flist.csv"
if [ ! -f "$fList" ]
then
  module load anaconda3
  conda activate /projects/above_gedi/users/pburns/envs/geospatial

  echo "----- Starting Script 0: Find and filter GEDI granules -----"
  srun --mem=5G --time=02:00:00 Rscript /home/pb463/scripts/repos/gedi_gridding/R/gedi_l2-l4a_s00_find_filter.R "$prod" "$region_name" "$region_geom" "$s_date" "$e_date" "$min_doy" "$max_doy" "$fList" "$l2a_dir" "$l2b_dir" "$l4a_dir"
  echo
fi


# ----- SLURM ARRAY -----
# Count files for SLURM Array
numFiles=$(wc -l "$fList" | cut -d " " -f1)

# Create the slurm array job script in the slurm folder
cd "${slurmDir}scripts"/
script_s01="sbatch_s01_${date_time}.sh"

# Specify the basename for the output log file
slurmLog_s01="${slurmDir}logs/s01_${date_time}_%A_%a.out"

cat > "$script_s01" <<EOT
#!/bin/bash

# - SLURM -
#SBATCH --job-name=gedil2l4_s01
#SBATCH --output="$slurmLog_s01"
#SBATCH --time="$timeReq_s1"
#SBATCH --mem="$memoryReq_s1"
#SBATCH --array=1-50000%"$arrayLim"
# NOTE: cant run more than 50000 jobs at once on our HPC, so need to chunk jobs manually

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

# Get the file names associated with each product
inFile_2a=\$(sed \$SLURM_ARRAY_TASK_ID'q;d' "$fList" | cut -d "," -f1)
inFile_2b=\$(sed \$SLURM_ARRAY_TASK_ID'q;d' "$fList" | cut -d "," -f2)
inFile_4a=\$(sed \$SLURM_ARRAY_TASK_ID'q;d' "$fList" | cut -d "," -f3)

# Run the script to extract, filter, join, and clip the GEDI products
echo "Starting RScript gedi_l2-l4a_s01_FJC.R"
echo
Rscript /home/pb463/scripts/repos/gedi_gridding/R/gedi_l2-l4a_s01_FJC.R \$inFile_2a \$inFile_2b \$inFile_4a "$region_name" "$region_geom" "$chunk_dim" "$s01Dir"

# - ENDING -
date_time=`date +%Y%m%d_%H%M%S`
echo "The ending date_time: " \$date_time
duration=\$SECONDS
echo "\$((\$duration / 60)) minutes and \$((\$duration % 60)) seconds elapsed."

EOT

# Run the slurm array job script
echo "----- Starting Script 1: Extract, filter, join, and chunk GEDI granules using SLURM arrays -----"
jobid_s01=$(sbatch --parsable "$script_s01")
echo "The SLURM job id is: ${jobid_s01}"


# Run jobstats after s01 is finished
slurmLog_s01js="${slurmDir}jobstats/s01js_${date_time}_${jobid_s01}.out"
cd "${slurmDir}scripts"/
script_s01js="sbatch_s01js_${date_time}.sh"

cat > "$script_s01js" <<EOT
#!/bin/bash

# - SLURM -
#SBATCH --job-name=gedil2l4_s01js
#SBATCH --output="$slurmLog_s01js"
#SBATCH --time=00:30:00
#SBATCH --mem=5G

srun sleep 5
srun /home/pb463/scripts/get_slurm_array_errors.sh "$jobid_s01"

EOT

jobid_s01js=$(sbatch --dependency=afterany:"$jobid_s01" "$script_s01js")
echo "jobstats will be run upon completion ${jobid_s01js}"
