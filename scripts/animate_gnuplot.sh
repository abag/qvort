#!/bin/bash
number_runs=$1
prefix=./data/var
number_mask=00
from=1
to=$number_runs
step=1
bsize=0.05
xmin=-$bsize
xmax=$bsize
ymin=-$bsize
ymax=$bsize
zmin=-$bsize
zmax=$bsize
for i in `seq $from $step $to`
do
  if [ $i -lt 10 ]
    then
    gnuplot - <<EOF
    set terminal png
    set output "${prefix}00$i.png"
    set xrange [${xmin}:${xmax}]
    set yrange [${ymin}:${ymax}]
    set zrange [${zmin}:${zmax}]
    sp "${prefix}000$i.log" u 1:2:3 w p
    set output
EOF
elif [ $i -lt 100 ]
    then
    gnuplot - <<EOF
    set terminal png
    set output "${prefix}0$i.png"
    set xrange [${xmin}:${xmax}]
    set yrange [${ymin}:${ymax}]
    set zrange [${zmin}:${zmax}]
    sp "${prefix}00$i.log" u 1:2:3 w p
    set output
EOF
elif [ $i -lt 1000 ]
    then
    gnuplot - <<EOF
    set terminal png
    set output "${prefix}$i.png"
    set xrange [${xmin}:${xmax}]
    set yrange [${ymin}:${ymax}]
    set zrange [${zmin}:${zmax}]
    sp "${prefix}0$i.log" u 1:2:3 w p
    set output
EOF
fi
done
animate data/*.png
