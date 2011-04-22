gnuplot << TOEND

#Some common options
set term dumb

plot 'data/ts.log' u 2:6 w l lt 1 lc 1 lw 2 t ''

TOEND
