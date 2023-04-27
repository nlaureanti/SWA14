#!/bin/bash


LS=$( ls *.nc.zip )

for l in $LS ;do
	echo $(echo $l | sed "s#.zip##g")
	#unzip $l $(echo $l | sed -i "s#.zip##g")
	mv data.nc $(echo $l | sed -i "s#.zip##g")
done

