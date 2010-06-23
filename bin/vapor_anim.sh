#!/bin/bash
  . /opt/vapor/bin/vapor-setup.sh
  filenum=$1
  nmesh=32
  rm -r u2*
  #vdfcreate -dimension $nmesh x$nmesh x$nmesh -numts $filenum -varnames u2 u2.vdf
  vdfcreate -dimension 64x64x64 -numts 100 -varnames u2 u2.vdf
  for i in `seq 1 $filenum`; do
    if [[ $i -le 9 ]] ; then
      file=data/vap_mesh00$i.dat
      echo $file
    else 
      file=data/vap_mesh0$i.dat
    fi
    raw2vdf -ts $i -varname u2 u2.vdf $file
  done
  vaporgui ./u2.vdf
