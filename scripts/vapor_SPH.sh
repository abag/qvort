#!/bin/bash
  . /opt/vapor/bin/vapor-setup.sh
  filenum=$1
  rm -rf ./SPH_rho*
  vdfcreate -dimension 64x64x64 -numts 100 -varnames rho SPH_rho.vdf
  for i in `seq 1 $filenum`; do
    if [[ $i -le 9 ]] ; then
      file1=data/SPH_vapor00$i.dat
    else 
      file1=data/SPH_vapor0$i.dat
    fi
    raw2vdf -ts $i -varname rho SPH_rho.vdf -dbl $file1
  done
  vaporgui 
