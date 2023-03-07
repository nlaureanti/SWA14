#!/bin/bash

dir=$(pwd)
exp='AScoast'

for v in 'temp' 'salt' 'speed' 'SSH' ; do
montage -tile x1 -mode concatenate -trim ${dir}/mom_${exp}_${v}_{1..5}.png ${dir}/MOM_${exp}_exp_${v}_dy.png
montage -tile x1 -mode concatenate -trim ${dir}/mom_${exp}_${v}_{1,15,30,58}.png ${dir}/MOM_${exp}_exp_${v}_pnt.png
 convert -delay 20 -trim +repage ${dir}/mom_${exp}_${v}_{1..50}.png ${dir}/MOM_${exp}_exp_${v}_movie.gif
done
