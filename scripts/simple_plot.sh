unset DISPLAY
matlab > matlab.out 2>&1 << EOF
octave_plot($1)
print -deps 'data/var$1.ps'
exit
EOF
