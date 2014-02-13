reset
set terminal epslatex standalone color colourtext font 'default, 14'
unset label
unset key
set origin 0.0,0.0
set size 1.0,1.0
set style line 1 lt 1 lc 7 lw 5
set output "figure.tex"
set xlabel '$t$'
set ylabel '$\langle \kappa \rangle$'
plot "<paste ../data/ts.log ../data/curvature.log " u 2:12 w l ls 1
set output
set term x11
