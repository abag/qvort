F90=gfortran
#FFLAGS=-O3 # Use 3rd level of optimisation
#FFLAGS=-O2 -Wall -fdefault-real-8           # DP/ 3rd level of optimisation pedantic
FFLAGS=-O2 -fdefault-real-8 -fopenmp    # DP/ 3rd level of optimisation
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
run.x: run.f90 cdata.o derivatives.o forcing.o reconnection.o initial_line.o initial_loop.o initial.o output.o periodic.o general.o statistics.o normal_fluid.o timestep.o line.o diagnostic.o quasip.o particles.o tree.o ksmodel.o mirror.o stiff_solver.o hamiltonian.o smoothing.o matrix.o topology.o inject.o kernel_density.o hyperviscous.o Makefile 
	$(F90) $(FFLAGS) -o ./run.x run.f90 cdata.o derivatives.o forcing.o normal_fluid.o reconnection.o initial.o initial_loop.o initial_line.o output.o periodic.o general.o statistics.o timestep.o line.o diagnostic.o quasip.o particles.o tree.o ksmodel.o mirror.o stiff_solver.o hamiltonian.o smoothing.o matrix.o topology.o inject.o kernel_density.o hyperviscous.o
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
statistics.o: statistics.f90 cdata.o
#
general.o: general.f90 derivatives.o statistics.o cdata.o
#
derivatives.o: derivatives.f90 cdata.o
#
normal_fluid.o: normal_fluid.f90 general.o cdata.o ksmodel.o
#
hyperviscous.o: hyperviscous.f90 general.o cdata.o
#
forcing.o: forcing.f90 general.o cdata.o
#
initial.o: initial.f90 initial_line.o initial_loop.o stiff_solver.o normal_fluid.o forcing.o cdata.o quasip.o periodic.o particles.o smoothing.o inject.o output.o hyperviscous.o
#
initial_loop.o: initial_loop.f90 cdata.o general.o periodic.o
#
initial_line.o: initial_line.f90 cdata.o general.o periodic.o
#
timestep.o: timestep.f90 normal_fluid.o forcing.o general.o cdata.o tree.o mirror.o matrix.o hyperviscous.o
#
line.o: line.f90 reconnection.o periodic.o general.o cdata.o matrix.o
#
reconnection.o: reconnection.f90 general.o cdata.o matrix.o
#
inject.o: inject.f90 general.o cdata.o
#
matrix.o: matrix.f90
#
diagnostic.o: diagnostic.f90 topology.o tree.o general.o cdata.o kernel_density.o
#
tree.o: tree.f90 general.o cdata.o
#
topology.o: topology.f90 general.o cdata.o
#
smoothing.o: smoothing.f90 tree.o general.o cdata.o
#
quasip.o: quasip.f90 hamiltonian.o general.o cdata.o timestep.o normal_fluid.o output.o tree.o stiff_solver.o
#
particles.o: particles.f90 general.o cdata.o timestep.o normal_fluid.o output.o tree.o
#
ksmodel.o: ksmodel.f90 general.o cdata.o
#
mirror.o: mirror.f90 general.o cdata.o
#
kernel_density.o: kernel_density.f90
#
stiff_solver.o: stiff_solver.f90 hamiltonian.o general.o cdata.o
#
hamiltonian.o: hamiltonian.f90 tree.o timestep.o output.o general.o cdata.o
#
