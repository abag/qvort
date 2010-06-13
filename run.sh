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
        -d|--double
                    Recompile the code using double precision.

        -q|--quiet
                    Use nohup to run the code in the background. This means 
                    the job will continue to run even if the terminal is closed.
                    All output is directed to a log file 'out.log'.
        -f|--force
                    If the script detects that data directory does not exist
                    then running with the -f option forces the creation of 
                    a directory, if on network space this can be slow...
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
options=`getopt -o compile,double,quiet,force,help -n run.sh -- "$@"`

# If no options, show the help
#if [ $# == 0 ]; then
#  usage
#  exit 1
#fi

eval set -- "$options"
COMPILE=0
DOUBLE=0
QUIET=0
RESTART=0
FORCE=0
# Set parameters depending on options
while true
do
  case "$1" in
    -c|--compile) COMPILE=1; shift;;
    -d|--double) DOUBLE=1; shift;;
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

echo Going to run...

if [ -d ./data ]; then
  echo ./data exists, proceeding...
else
  if [ $FORCE -eq 0 ]; then
    echo WARNING: data directory does not exist - aborting.
    echo Use -f to force the creation of a directory, or create it yourself.
    exit 1
  else
    mkdir data
  fi
fi
if [ $RESTART -eq 1 ]; then
    echo code will restart if possible
  else
  if [ -f ./data/var.dat ]; then
    echo deleting varfile
    rm ./data/var.dat 
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

if [ $QUIET -eq 1 ]; then
  # Run quietly
  run_quiet
  exit 0
fi

# Run the job.
run

