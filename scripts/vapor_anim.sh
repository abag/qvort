#!/bin/bash
  . /opt/vapor/bin/vapor-setup.sh
  filenum=$1
  nmesh=32
  rm -r u2*
  #vdfcreate -dimension $nmesh x$nmesh x$nmesh -numts $filenum -varnames u2 u2.vdf
  #vdfcreate -dimension 64x64x64 -numts 100 -varnames u2 u2.vdf
  vdfcreate -dimension 128x128x128 -numts 100 -varnames unorm unorm.vdf
  vdfcreate -dimension 128x128x128 -numts 100 -varnames usup usup.vdf
  for i in `seq 1 $filenum`; do
    if [[ $i -le 9 ]] ; then
      file1=data/vap_norm00$i.dat
      file2=data/vap_sup00$i.dat
    else 
      file1=data/vap_norm0$i.dat
      file2=data/vap_sup0$i.dat
    fi
    raw2vdf -ts $i -varname unorm unorm.vdf -dbl $file1
    raw2vdf -ts $i -varname usup usup.vdf -dbl $file2
  done
  vaporgui 
