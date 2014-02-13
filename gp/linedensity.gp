reset
set terminal epslatex standalone color colourtext font 'default, 14'
unset label
unset key
set origin 0.0,0.0
set size 1.0,1.0
set style line 1 lt 1 lc 7 lw 5
set output "figure.tex"
set xlabel '$t$'
set ylabel '$L$'
D = "`sed -n '2p' ../data/dims.log | awk '{print $1}'`"
p '../data/ts.log' u ($2):($6/D**3) w l ls 1
set output
set term x11
