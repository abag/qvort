reset
set terminal epslatex standalone color colourtext font 'default, 12'
unset label
unset key
set origin 0.0,0.0
set size 1.0,1.0
set style line 1 lt 1 lc 7 lw 5
set output "figure.tex"
set xlabel '$t$'
set ylabel '$\langle u_\mathrm{BS} \rangle / \langle u_\mathrm{LIA} \rangle$'
plot '../data/local_v_nonlocal.log' u 2:($3) w l ls 1
set output
set term x11
