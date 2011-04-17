#!/bin/bash
usage() {
cat <<EOF
qvort_nr.sh
---------

DESCRIPTION:

Clone the qvort code to another directory

USAGE:

        run.sh [OPTIONS] new_directory

OPTIONS:        
        -h
                    Show usage information.
EOF
}
# Get command line options
options=`getopt -o compile,protect,quiet,force,help -n run.sh -- "$@"`
# If no options, show the help
if [ $# == 0 ]; then
  echo "No destination directory set, showing help"
  usage
  exit 1
fi
eval set -- "$options"
while true
do
  case "$1" in
    -h|--help) usage; exit 1;;
    --) shift ; break ;;
    *) echo "Invalid flag" ; exit 1 ;;
  esac
done
#THE MOST IMPORTANT LINE IN THIS SCRIPT MAKE SURE IT IS SET CORRECTLY
PRJ_HOME=${HOME}/Work/code/qvort
#i will check it at least exists
if [ ! -d "$PRJ_HOME" ]; then
  echo "arghhhh, qvort home directory does not exist!!"
  exit 1
fi
#now I will check if it is incorrectly assigned
if [ ! -f "$PRJ_HOME/Doxyfile" ] ; then
  echo "the directory you have set as qvort's home exist but is incorrect..."
  exit 1
fi
if [ ! -f "$PRJ_HOME/run.sh" ] ; then
  echo "the directory you have set as qvort's home exist but is incorrect..."
  exit 1
fi
if [ ! -f "$PRJ_HOME/README" ] ; then
  echo "the directory you have set as qvort's home exist but is incorrect..."
  exit 1
fi
if [ ! -f "$PRJ_HOME/run.in" ] ; then
  echo "the directory you have set as qvort's home exist but is incorrect..."
  exit 1
fi
if [ ! -f "$PRJ_HOME/startup.m" ] ; then
  echo "the directory you have set as qvort's home exist but is incorrect..."
  exit 1
fi
#where are we?
startdir=`pwd`
rundir=
#check the rundir doesn't exists already
if [ "$1" ]; then
  rundir="$1"
  if [ -d "$rundir" ]; then
    echo "Rundir '$rundir' already exists!!"
    exit 1
  fi
  mkdir "$rundir"
  cd "$rundir"
fi
#make src if needed
[ -d src ] || mkdir src
cd src
#link all src files here
ln -s ${PRJ_HOME}/src/Makefile ${PRJ_HOME}/src/*.f90 .
#move out of src 
cd ..
#makefile
ln -s ${PRJ_HOME}/Makefile Makefile 
#startup script
ln -s ${PRJ_HOME}/run.sh .
#run.in must cp not ln, this is unique to this run
cp ${PRJ_HOME}/run.in .
#make the data directory
if [ ! -d "data" ]; then
  mkdir data
fi
#finally create a new startup.m file which contains link to PRJ_HOME
if [ ! -f "startup.m" ]; then
  echo "%to plot a single snapshot run vortex_plot(filenumber)" > startup.m
  echo "%e.g. vortex_plot(1)" >> startup.m
  echo "%to get time series information run ts" >>startup.m
  echo "disp('adding all qvort paths')" >>startup.m
  echo "addpath(genpath('$PRJ_HOME'))" >>startup.m
  echo "disp('type help startup to get simple plotting commands for qvort')" >>startup.m
  echo "successfully cloned qvort to $rundir"
fi

