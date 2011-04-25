#!/bin/bash
unset DISPLAY
matlab > matlab.out 2>&1 << EOF
vortex_anim($1,$2,2)
exit
EOF
