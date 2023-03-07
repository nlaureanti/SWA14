#!/bin/bash

dir=$(pwd)
exp='SWA002'

for v in 'temp' 'salt' 'speed' 'SSH' ; do
montage -tile x1 -mode concatenate -trim ${dir}/${v}/mom_${exp}_${v}_sup_{20011103..20011108}0000.png ${dir}/MOM_${exp}_exp_${v}_dy.png
montage -tile x1 -mode concatenate -trim ${dir}/${v}/mom_${exp}_${v}_sup_{20011130,20011203,20011215,20011228}0000.png ${dir}/MOM_${exp}_exp_${v}_pnt.png
montage -tile x1 -mode concatenate -trim ${dir}/${v}/mom_${exp}_${v}_sup_{20011203,20011208,20011214,20011219}0000.png ${dir}/MOM_${exp}_exp_${v}_SACZinac.png
# convert -delay 20 -trim +repage ${dir}/mom_${exp}_${v}_{1..50}.png ${dir}/MOM_${exp}_exp_${v}_movie.gif
done
