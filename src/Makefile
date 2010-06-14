BUILD=$1
F90=gfortran
#FFLAGS=-fdefault-real-8            # Double precision
FFLAGS=-O3           # Use 3rd level of optimisation
DFLAGS=-O3 -fdefault-real-8           # DP/ 3rd level of optimisation
   #
   #  -pg   =  profiler (gprof)
   #  -g    = debugger  (gdb)
   #

## Generic rules:
.SUFFIXES:          # get rid of that annoying Modula rule
.SUFFIXES: .f .f90 .c .o .obj .x .exe .inc .h  

#
# Rules (how to make a Blergh from a Blah
#
ifeq "$@" "double"
.f90.o:
	    $(F90) $(DFLAGS) -o $*.o -c $*.f90
else
.f90.o:
	    $(F90) $(FFLAGS) -o $*.o -c $*.f90
endif

#
# Special targets
#
default: run.x
#
# Targets
#
run.x: run.f90 cdata.o initial.o output.o periodic.o general.o timestep.o line.o diagnostic.o quasip.o Makefile 
	$(F90) $(FFLAGS) -o ./run.x run.f90 cdata.o initial.o output.o periodic.o general.o timestep.o line.o diagnostic.o quasip.o
double: run.f90 cdata.o initial.o output.o periodic.o general.o timestep.o line.o diagnostic.o quasip.o Makefile
	$(F90) $(DFLAGS) -o ./run.x run.f90 cdata.o initial.o output.o periodic.o general.o timestep.o line.o diagnostic.o quasip.o
#
#
clean: 
	rm -f *.o *.x *.exe *.obj *.mod
cleann: clean
	rm -f *~
pristine: cleann
	rm -f *.dat *.log *.png
#
# Dependencies
#
cdata.o: cdata.f90 ../run.in
# 
initial.o: initial.f90 cdata.o
#
output.o: output.f90 cdata.o
#
periodic.o: periodic.f90 cdata.o
#
general.o: general.f90 cdata.o
#
timestep.o: timestep.f90 general.o cdata.o
#
line.o: line.f90 general.o cdata.o
#
diagnostic.o: diagnostic.f90 general.o cdata.o
#
quasip.o: quasip.f90 general.o cdata.o
