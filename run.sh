#!/bin/bash
#$Id: run.sh 111 2010-05-16 15:02:54Z abag $
usage() {
cat <<EOF
run.sh
---------

DESCRIPTION:

Run the qvort code

USAGE:

        run.sh [OPTIONS]

OPTIONS:        

        -c|--compile
                    Recompile the code - useful if run.in has been changed.
        -p|--protect
                    If --force is used the protect flag protects the random
                    seed from being removed
        -q|--quiet
                    Use nohup to run the code in the background. This means 
                    the job will continue to run even if the terminal is closed.
                    All output is directed to a log file 'out.log'.
        -f|--force
                    If the script detects that data directory is not empty
                    then running with the -f option will empty the data
                    directory
        -r|--restart
                    Set this option to restart the code from the last store Var
                    file
        -h|--help
                    Show usage information.
EOF
}
run() {
  echo Running job...
  /usr/bin/time -p $EXE
}
run_quiet() {
  echo Running job quietly...
  nohup /usr/bin/time -p $EXE &> $LOGFILE &
}

# Get command line options
options=`getopt -o compile,protect,quiet,force,help -n run.sh -- "$@"`

# If no options, show the help
#if [ $# == 0 ]; then
#  usage
#  exit 1
#fi

eval set -- "$options"
COMPILE=0
PROTECT=0
QUIET=0
RESTART=0
FORCE=0
# Set parameters depending on options
while true
do
  case "$1" in
    -c|--compile) COMPILE=1; shift;;
    -p|--protect) PROTECT=1; shift;;
    -q|--quite) QUIET=1; shift;;
    -f|--force) FORCE=1; shift;;
    -r|--restart) RESTART=1; shift;;
    -h|--help) usage; exit 1;;
    --) shift ; break ;;
    *) echo "Invalid flag" ; exit 1 ;;
  esac
done

#if [ $# -eq 0 ];
#then
#  usage
#  exit 1
#fi

EXE=./src/run.x
LOGFILE=out.log


if [ -d ./data ]; then
  echo ./data exists, proceeding...
else
  echo "WARNING: data directory does not exist - aborting."
  echo "create the directory using mkdir data."
  exit 1
fi
if [ $FORCE -eq 1 ]; then
  if [ $PROTECT -eq 1 ]; then
    echo "emptying data but protecting random seed"
    if [ -e ./data/seed.dat ]; then
      mv data/seed.dat .
      rm -f data/*
      mv seed.dat ./data
    else
      echo "there is no seed to protect!"
    fi
  else
    echo "emptying data"
    rm -f data/*
  fi
  if [ -e ./STOP ]; then
    echo "removing STOP file"
    rm ./STOP
  fi
fi
if [ -e ./STOP ]; then
  echo "encountered STOP file: aborting run"
  echo "please delete the STOP file before running again"
  exit 1
fi
if [ $RESTART -eq 1 ]; then
  echo code will restart if possible
else
  if [ -f ./data/var.dat ]; then
    echo emptying data directory
    rm ./data/* 
  fi
fi
if [ $COMPILE -eq 1 ]; then
  # Recompile
  make clean
  if [ $DOUBLE -eq 1 ]; then
    echo Recompiling code in double precision...
    make double
  else
    echo Recompiling code...
    make
  fi
fi

#check that executable file exists
if [ ! -e $EXE ] ; then
  echo "WARNING: exectuable does not exist"
  echo "please compile the code and re-run"
  exit 1
fi

echo -e "\033[1mGoing to run...\033[0m"

if [ $QUIET -eq 1 ]; then
  # Run quietly
  run_quiet
  exit 0
fi

# Run the job.
run

