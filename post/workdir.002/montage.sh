#!/bin/bash

dir=$(pwd)
exp="SWA002${1}"

vvars=$( ls $dir )


#IC jan
for v in $vvars ; do
	ls ${dir}/${v}/mom_${exp}_${v}_* >/dev/null | exit 0
montage -tile x1 -mode concatenate -trim ${dir}/${v}/mom_${exp}_${v}_sup_{20010103..20010108}0000.png ${dir}/MOM_${exp}_exp_${v}_dy.png
montage -tile x1 -mode concatenate -trim ${dir}/${v}/mom_${exp}_${v}_sup_{20010103,20010108,20010113,20010118}0000.png ${dir}/MOM_${exp}_exp_${v}_5dy.png
montage -tile x1 -mode concatenate -trim ${dir}/${v}/mom_${exp}_${v}_sup_{20010103,20010118,20010202,20010217}0000.png ${dir}/MOM_${exp}_exp_${v}_pnt.png
#montage -tile x1 -mode concatenate -trim ${dir}/${v}/mom_${exp}_${v}_sup_{20010203,20010208,20011214,20011219}0000.png ${dir}/MOM_${exp}_exp_${v}_SACZinac.png
# convert -delay 20 -trim +repage ${dir}/mom_${exp}_${v}_{1..50}.png ${dir}/MOM_${exp}_exp_${v}_movie.gif
done


exit

#IC nov:::
for v in $vvars ; do
montage -tile x1 -mode concatenate -trim ${dir}/${v}/mom_${exp}_${v}_sup_{20011103..20011108}0000.png ${dir}/MOM_${exp}_exp_${v}_dy.png
montage -tile x1 -mode concatenate -trim ${dir}/${v}/mom_${exp}_${v}_sup_{20011130,20011203,20011215,20011228}0000.png ${dir}/MOM_${exp}_exp_${v}_pnt.png
montage -tile x1 -mode concatenate -trim ${dir}/${v}/mom_${exp}_${v}_sup_{20011203,20011208,20011214,20011219}0000.png ${dir}/MOM_${exp}_exp_${v}_SACZinac.png
# convert -delay 20 -trim +repage ${dir}/mom_${exp}_${v}_{1..50}.png ${dir}/MOM_${exp}_exp_${v}_movie.gif
done
