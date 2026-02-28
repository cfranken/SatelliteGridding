#!/bin/bash i
# Directory where NetCDF files are located 
DIR="/net/fluo/data2/projects/TROPOMI_GASES/ch4/" 
# Loop over all NetCDF files in the directory 
for file in $DIR/*.nc; do 
	# Use ncdump to check the file 
	ncdump -h $file > /dev/null 2>&1 
	# If ncdump fails, then delete the file 
	if [ $? -ne 0 ]; then 
		echo "Deleting corrupt file: $file" 
		rm $file 
	fi 
done
