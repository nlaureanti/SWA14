#!/bin/bash

pwd=$(pwd)
ls=$(ls $pwd/*.nc) 

mkdir classic

for f in $ls ;do
        echo  classic/$(basename $f  )
	ncks -7 $f classic/$(basename $f  )
	
done
