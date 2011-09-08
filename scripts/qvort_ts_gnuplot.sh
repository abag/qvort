gnuplot << TOEND
set terminal postscript eps color enhanced
set output 'ts.eps'

#Some common options
set grid

#This begins the multiplot environment. This example uses a 2x2 regular layout
set multiplot layout 3,3 title 'Filament time series plots'

set xtics font "Arial,6"
set ytics font "Arial,6"
set size ratio 0.5

set title 'point count'
plot './data/ts.log' u 3 w l lt 1 lc 1 lw 2 t ''

set title 'recon count'
plot './data/ts.log' u 4 w l lt 1 lc 2 lw 2 t ''

set title 'avg sep'
plot './data/ts.log' u 5 w l lt 1 lc 3 lw 2 t ''

set title 'filament length'
plot './data/ts.log' u 6 w l lt 1 lc 4 lw 2 t ''

set title 'max u'
plot './data/ts.log' u 7 w l lt 1 lc 5 lw 2 t ''

set title 'max du'
plot './data/ts.log' u 8 w l lt 1 lc 6 lw 2 t ''

set title '# evaluations per point'
plot './data/ts.log' u 9 w l lt 1 lc 7 lw 2 t ''

set title 'mean curv'
plot './data/curvature.log' u 1 w l lt 1 lc 8 lw 2 t ''

set title 'max curv'
plot './data/curvature.log' u 3 w l lt 1 lc 9 lw 2 t ''

#close the multiplot environment before exit
unset multiplot

TOEND
