#!/bin/bash

#######################################################
# Search for Data
#https://oceandata.sci.gsfc.nasa.gov/api/file_search/

#wget -q --post-data="results_as_file=1&sensor_id=6&dtid=1123&sdate=1997-09-04 00:00:00&edate=2010-12-11 00:00:00&subType=1&addurl=1&prod_id=chlor_a&resolution_id=9km&period=DAY" -O - https://oceandata.sci.gsfc.nasa.gov/api/file_search

#Before download, use app key
#https://oceancolor.gsfc.nasa.gov/data/download_methods/#appkey
######################################################

file_to_download=$1

if [[ ${#1} == 0 ]] ;then
read -p "insert file_to_download    " file_to_download
fi

string="?appkey=d9e257104213052988ddf3204cc59ad6febf28e6"
if [ -f $file_to_download ] ; then
	for line in $(cat $file_to_download) ; do
		wget ${line}${string} 
	done
fi

cdo='/home/nicole/workdir/envs/cdo_env/bin/cdo'
mkdir tmp
for file in $(ls SEASTAR_SEAWIFS_GAC.*) ; do
	datestr=$(echo ${file}| cut -d "." -f 2)
	echo ${datestr:0:4}-${datestr:4:2}-${datestr:6:8}
	date_to_cdo=${datestr:0:4}-${datestr:4:2}-${datestr:6:8}
	new_f=$(echo ${file}| cut -d "." -f 1)_${datestr:0:4}-${datestr:4:2}-${datestr:6:8}.nc
	echo $new_f
#	$cdo -setdate,${datestr:0:4}-${datestr:4:2}-${datestr:6:8} $file $new_f
	$cdo  -settaxis,$date_to_cdo,00:00:00,1day $file $new_f
	mv $file tmp
	
done
date
rm -f ./SEAWIFSDAILY9KM.nc
$cdo mergetime SEASTAR_SEAWIFS_GAC*.nc ./SEAWIFSDAILY9KM.nc
date

echo DONE!
