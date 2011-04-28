#!/bin/bash
number_runs=$1
prefix=./data/var
number_mask=00
from=1
to=$number_runs
step=1
#read in the box dimensions
i=0
while read LINE
do
    i=$(($i+1))
    if [ "$i" -eq 2 ] ; then
    bsize=$LINE
    echo $bsize
    fi
done < ./data/dims.log
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
    set xrange [${xmin}/2:${xmax}/2]
    set yrange [${ymin}/2:${ymax}/2]
    set zrange [${zmin}/2:${zmax}/2]
    sp "${prefix}000$i.log" u 1:2:3 w p
    set output
EOF
elif [ $i -lt 100 ]
    then
    gnuplot - <<EOF
    set terminal png
    set output "${prefix}0$i.png"
    set xrange [${xmin}/2:${xmax}/2]
    set yrange [${ymin}/2:${ymax}/2]
    set zrange [${zmin}/2:${zmax}/2]
    sp "${prefix}00$i.log" u 1:2:3 w p
    set output
EOF
elif [ $i -lt 1000 ]
    then
    gnuplot - <<EOF
    set terminal png
    set output "${prefix}$i.png"
    set xrange [${xmin}/2:${xmax}/2]
    set yrange [${ymin}/2:${ymax}/2]
    set zrange [${zmin}/2:${zmax}/2]
    sp "${prefix}0$i.log" u 1:2:3 w p
    set output
EOF
fi
done
animate data/*.png
