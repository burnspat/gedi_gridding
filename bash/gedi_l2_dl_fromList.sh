#!/bin/bash

## ----- ABOUT -----
# gedil2_dl_fromList.sh
# Author: Patrick Burns [pb463@nau.edu]
# About: Download all (or subsets of) GEDI data from NASA Datapool using a list of URLS
# Last update: 7 Sept.2022
# TODO:
#



# ----- INPUTS -----

# GEDI specifics
# download directory
dlDir=$1
#'/projects/geode_data/share/lidar/GEDI02_B.002/0_orig/'

# file list to download
wget_flist_i=$2
#$dlDir'g2b002_global_20220412thru20221221.txt'


# NASA Earthdata specifics
#read -p 'NASA Earthdata Username: ' uservar
#echo
#read -sp 'NASA Earthdata Password: ' passvar
#echo
uservar=$3
passvar=$4

# ----- PROCESSING -----

echo
echo "Starting to process request..."
echo
# Access date
accDate=$(date '+%Y%m%d')

# Make sure there are no carriage returns
wget_base=$(basename "$wget_flist_i" .txt)
wget_flist="${dlDir}${wget_base}"_ncr.txt
if [ ! -f "$wget_flist" ]
then
  sed 's/\r$//' "$wget_flist_i" > "$wget_flist"
  rm "$wget_flist_i"
fi

# Make processing directory structure
# Check if directories already exist, if not make them
if [ ! -d "$dlDir" ]
then
    mkdir -p "$dlDir"
    echo "Download directory DID NOT exist. Made directory $dlDir"
    echo
fi

alldataDir="$dlDir"alldata/
if [ ! -d "$alldataDir" ]
then
    mkdir -p "$alldataDir"
    mkdir -p "$dlDir"alldata/incomplete/
fi

indexDir="$dlDir"index/
if [ ! -d "$indexDir" ]
then
    mkdir -p "$indexDir"
fi

cd "$alldataDir"

# Output an index of downloaded files
base=$(basename "$wget_flist" | cut -d "." -f1)
outIndex="${indexDir}${base}"_dl_status.csv
if [ -f "$outIndex" ]
then
    # Remove old completed index list
    rm "$outIndex"

    # Make new completed index list
    touch "$outIndex"
else
    # Make new completed index list
    touch "$outIndex"
fi

# Start downloading h5 files.
START=1
numFiles=$(wc -l "$wget_flist" | cut -d " " -f1)

echo "$numFiles files to download. Now might be a good time to check your storage capacity"
echo

for (( i="$START"; i<=(( "$numFiles" + 1 )); i++ ))
do
    dlPath=$(sed "$i"'q;d' "$wget_flist")
    dlPath_x="$dlPath".xml
    dlBase=$(basename "$dlPath" .h5)
    dlFile="${alldataDir}${dlBase}".h5
    dlFile_x="$dlFile".xml

    echo "Working on: $dlPath"
    echo "dlPath_x: $dlPath_x"
    echo "dlBase: $dlBase"
    echo "dlFile: $dlFile"
    echo "dlFile_x: $dlFile_x"

    # Check if file already exists in the specified download directory
    if [ ! -f "$dlFile" ]
    then
        echo "H5 file does not exist in specified download directory."
        # If file doesn't exist then download it
        echo "Downloading H5: "
        # Download first to a staging folder "incomplete/"
        wget -q --show-progress --no-parent --no-directories --continue -P "$alldataDir"incomplete/ --http-user=$uservar --http-password="$passvar" "$dlPath"

        # Once download has finished, move it
        mv "$alldataDir"incomplete/"$dlBase".h5 "$dlFile"

        # Also download associated xml
        # Remove existing file if necessary
        if [ -f "$dlFile_x" ]
        then
          rm "$dlFile_x"
        fi

        echo "Downloading XML:"
        wget -q --show-progress --no-parent --no-directories --continue -P "$alldataDir"incomplete/ --http-user="$uservar" --http-password="$passvar" "$dlPath_x"

        # Once xml download has finished, move it
        mv "$alldataDir"incomplete/"$dlBase".h5.xml "$dlFile_x"

        echo "Verifying checksums:"
        # Get the md5sum from the xml
        md5_x=$(grep 'Checksum>' "$dlFile_x" | cut -d ">" -f2 | cut -d "<" -f1)
        echo "XML checksum: $md5_x"

        # Get the md5sum from the file
        md5_f=$(md5sum "$dlFile" | cut -d " " -f1)
        echo " H5 checksum: $md5_f"

        # Check to see if md5 from xml matches md5 from file
        if [ "$md5_x" = "$md5_f" ]
        then
          echo "Checksums match!"
          chk1=1

          # Save status to .csv
          out_line="$dlPath",new_download,try1,"$chk1"

          # Add this file name to a completed list
          echo "$out_line" >> "$outIndex"
        else
          echo "Checksums do NOT match!"
          echo "Removing original files."

          # Remove original files since checksums don't match
          rm "$dlFile"
          rm "$dlFile_x"

          echo "Re-downloading H5:"
          # Download first to a staging folder "incomplete/"
          wget -q --show-progress --no-parent --no-directories --continue -P "$alldataDir"incomplete/ --http-user=$uservar --http-password="$passvar" "$dlPath"

          # Once download has finished, move it
          mv "$alldataDir"incomplete/"$dlBase".h5 "$dlFile"

          # Also download associated xml
          echo "Re-downloading XML:"
          wget -q --show-progress --no-parent --no-directories --continue -P "$alldataDir"incomplete/ --http-user="$uservar" --http-password="$passvar" "$dlPath_x"

          # Once xml download has finished, move it
          mv "$alldataDir"incomplete/"$dlBase".h5.xml "$dlFile_x"

          echo "Verifying checksums (round 2):"
          # Get the md5sum from the xml
          md5_x=$(grep 'Checksum>' "$dlFile_x" | cut -d ">" -f2 | cut -d "<" -f1)
          echo "XML checksum: $md5_x"

          # Get the md5sum from the file
          md5_f=$(md5sum "$dlFile" | cut -d " " -f1)
          echo " H5 checksum: $md5_f"

          # Check to see if md5 from xml matches md5 from file
          if [ "$md5_x" = "$md5_f" ]
          then
            echo "Checksums match!"
            chk2=1
          else
            echo "Checksums do NOT match after 2 tries!"
            chk2=0
          fi

          # Save status to .csv
          out_line="$dlPath",new_download,try2,"$chk2"

          # Add this file name to a completed list
          echo "$out_line" >> "$outIndex"

        fi

        echo "Done with $dlBase".h5

        # Calculate and report percent progress
        percDone=$(awk "BEGIN {print 100*$i/$numFiles}")
        echo "Overall progress is $percDone%"
        echo
    else
        echo "H5 File already exists in specified download directory."

        # Remove existing .xml file if necessary
        if [ -f "$dlFile_x" ]
        then
          rm "$dlFile_x"
        fi

        # Re-download the .xml file
        echo "Downloading XML:"
        wget -q --show-progress --no-parent --no-directories --continue -P "$alldataDir"incomplete/ --http-user="$uservar" --http-password="$passvar" "$dlPath_x"

        # Once xml download has finished, move it
        mv "$alldataDir"incomplete/"$dlBase".h5.xml "$dlFile_x"

        echo "Verifying checksums:"
        # Get the md5sum from the xml
        md5_x=$(grep 'Checksum>' "$dlFile_x" | cut -d ">" -f2 | cut -d "<" -f1)
        echo "XML checksum: $md5_x"

        # Get the md5sum from the file
        md5_f=$(md5sum "$dlFile" | cut -d " " -f1)
        echo " H5 checksum: $md5_f"

        # Check to see if md5 from xml matches md5 from file
        if [ "$md5_x" = "$md5_f" ]
        then
          echo "Checksums match!"
          chk1=1
          # Save status to .csv
          out_line="$dlPath",old_download,try1,"$chk1"

          # Add this file name to a completed list
          echo "$out_line" >> "$outIndex"

        else
          echo "Checksums do NOT match!"
          echo "Removing original files."

          # Remove original files since checksums don't match
          rm "$dlFile"
          rm "$dlFile_x"

          # Download first to a staging folder "incomplete/"
          echo "Re-downloading H5:"
          wget -q --show-progress --no-parent --no-directories --continue -P "$alldataDir"incomplete/ --http-user="$uservar" --http-password="$passvar" "$dlPath"

          # Once download has finished, move it
          mv "$alldataDir"incomplete/"$dlBase".h5 "$dlFile"

          # Also download associated xml
          echo "Re-downloading XML:"
          wget -q --show-progress --no-parent --no-directories --continue -P "$alldataDir"incomplete/ --http-user="$uservar" --http-password="$passvar" "$dlPath_x"

          # Once xml download has finished, move it
          mv "$alldataDir"incomplete/"$dlBase".h5.xml "$dlFile_x"

          echo "Verifying checksums (round 2):"
          # Get the md5sum from the xml
          md5_x=$(grep 'Checksum>' "$dlFile_x" | cut -d ">" -f2 | cut -d "<" -f1)
          echo "XML checksum: $md5_x"

          # Get the md5sum from the file
          md5_f=$(md5sum $dlFile | cut -d " " -f1)
          echo " H5 checksum: $md5_f"

          # Check to see if md5 from xml matches md5 from file
          if [ "$md5_x" = "$md5_f" ]
          then
            echo "Checksums match!"
            chk2=1
          else
            echo "Checksums do NOT match after 2 tries!"
            chk2=0
          fi

          # Save status to .csv
          out_line="$dlPath",old_download,try2,"$chk2"

          # Add this file name to a completed list
          echo "$out_line" >> "$outIndex"
        fi

        echo "Done with $dlBase".h5

        # Calculate and report percent progress
        percDone=$(awk "BEGIN {print 100*$i/$numFiles}")
        echo "Overall progress is $percDone%"
        echo
    fi
done
