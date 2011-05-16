#!/bin/bash
usage() {
cat <<EOF
qvort_relink.sh
---------

DESCRIPTION:

Relink files in src after a git update

USAGE:

        qvort_relink [OPTIONS]

OPTIONS:        
        -h
                    Show usage information.
EOF
}
# Get command line options
options=`getopt -o help -n qvort_relink -- "$@"`
# If no options, show the help
if [ $# == 0 ]; then
  #echo "No destination directory set, showing help"
  #usage
  #exit 1
  rundir=`pwd`
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
PRJ_HOME=${HOME}/svnprojects/qvort
#i will check it at least exists
if [ ! -d "$PRJ_HOME" ]; then
  echo "arghhhh, qvort home directory does not exist!!"
  exit 1
fi
#now I will check if it is incorrectly assigned
if [ ! -f "$PRJ_HOME/Doxyfile" ] ; then
  echo "the directory you have set as qvort's home exists but is incorrect..."
  exit 1
fi
if [ ! -f "$PRJ_HOME/run.sh" ] ; then
  echo "the directory you have set as qvort's home exists but is incorrect..."
  exit 1
fi
if [ ! -f "$PRJ_HOME/README" ] ; then
  echo "the directory you have set as qvort's home exists but is incorrect..."
  exit 1
fi
if [ ! -f "$PRJ_HOME/run.in" ] ; then
  echo "the directory you have set as qvort's home exists but is incorrect..."
  exit 1
fi
if [ ! -f "$PRJ_HOME/startup.m" ] ; then
  echo "the directory you have set as qvort's home exists but is incorrect..."
  exit 1
fi
#where are we?
startdir=`pwd`
rundir=
#check the rundir exists
if [ "$1" ]; then
  rundir="$1"
  if [ ! -d "$rundir" ]; then
    echo "Directory '$rundir' doesnt exist!"
    exit 1
  fi
fi
cd src
#link all src files here
ln -sf ${PRJ_HOME}/src/Makefile ${PRJ_HOME}/src/*.f90 .
#move out of src 
cd ..
echo "successfully relinked files in src"

