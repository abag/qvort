F90=gfortran
#FFLAGS=-O3 # Use 3rd level of optimisation
FFLAGS=-O3 -fdefault-real-8           # DP/ 3rd level of optimisation
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
.f90.o:
	    $(F90) $(FFLAGS) -o $*.o -c $*.f90

#
# Special targets
#
default: run.x
#
# Targets
#
run.x: run.f90 cdata.o forcing.o initial.o output.o periodic.o general.o normal_fluid.o timestep.o line.o diagnostic.o quasip.o tree.o Makefile 
	$(F90) $(FFLAGS) -o ./run.x run.f90 cdata.o forcing.o normal_fluid.o initial.o output.o periodic.o general.o timestep.o line.o diagnostic.o quasip.o tree.o
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
#
output.o: output.f90 cdata.o tree.o
#
periodic.o: periodic.f90 cdata.o
#
general.o: general.f90 cdata.o
#
normal_fluid.o: normal_fluid.f90 general.o cdata.o
#
forcing.o: forcing.f90 general.o cdata.o
#
initial.o: initial.f90 normal_fluid.o forcing.o cdata.o quasip.o
#
timestep.o: timestep.f90 normal_fluid.o forcing.o general.o cdata.o tree.o
#
line.o: line.f90 general.o cdata.o
#
diagnostic.o: diagnostic.f90 general.o cdata.o
#
tree.o: tree.f90 general.o cdata.o
#
quasip.o: quasip.f90 general.o cdata.o timestep.o normal_fluid.o output.o
