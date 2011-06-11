#!/bin/bash
  . /opt/vapor/bin/vapor-setup.sh
  vdfcreate -dimension 64x64x64 -numts 1 -varnames var snap.vdf
  raw2vdf -ts 0 -varname var snap.vdf -dbl $1
  vaporgui 
