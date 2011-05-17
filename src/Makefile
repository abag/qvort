F90=gfortran
#FFLAGS=-O3 # Use 3rd level of optimisation
#FFLAGS=-O2 -Wall -fdefault-real-8           # DP/ 3rd level of optimisation pedantic
FFLAGS=-O2 -fdefault-real-8     # DP/ 3rd level of optimisation
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
run.x: run.f90 cdata.o forcing.o initial_cond.o initial.o output.o periodic.o general.o normal_fluid.o timestep.o line.o diagnostic.o quasip.o particles.o tree.o ksmodel.o mirror.o stiff_solver.o hamiltonian.o smoothing.o mag.o sph.o sph_tree.o sph_kernel.o matrix.o sph_interface.o topology.o inject.o Makefile 
	$(F90) $(FFLAGS) -o ./run.x run.f90 cdata.o forcing.o normal_fluid.o initial.o initial_cond.o output.o periodic.o general.o timestep.o line.o diagnostic.o quasip.o particles.o tree.o ksmodel.o mirror.o stiff_solver.o hamiltonian.o smoothing.o mag.o sph.o sph_tree.o sph_kernel.o matrix.o sph_interface.o topology.o inject.o 
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
normal_fluid.o: normal_fluid.f90 general.o cdata.o ksmodel.o
#
forcing.o: forcing.f90 general.o cdata.o
#
initial.o: initial.f90 initial_cond.o stiff_solver.o normal_fluid.o forcing.o cdata.o quasip.o periodic.o particles.o smoothing.o sph.o inject.o
#
initial_cond.o: initial_cond.f90 cdata.o general.o periodic.o
#
timestep.o: timestep.f90 normal_fluid.o forcing.o general.o cdata.o tree.o mirror.o mag.o matrix.o sph_interface.o
#
line.o: line.f90 general.o cdata.o
#
inject.o: inject.f90 general.o cdata.o
#
matrix.o: matrix.f90
#
diagnostic.o: diagnostic.f90 topology.o general.o cdata.o
#
tree.o: tree.f90 general.o cdata.o
#
topology.o: topology.f90 general.o cdata.o
#
smoothing.o: smoothing.f90 tree.o general.o cdata.o
#
mag.o: mag.f90 general.o cdata.o smoothing.o tree.o
#
quasip.o: quasip.f90 hamiltonian.o general.o cdata.o timestep.o normal_fluid.o output.o tree.o stiff_solver.o
#
particles.o: particles.f90 general.o cdata.o timestep.o normal_fluid.o output.o tree.o
#
ksmodel.o: ksmodel.f90 general.o cdata.o
#
mirror.o: mirror.f90 general.o cdata.o
#
stiff_solver.o: stiff_solver.f90 hamiltonian.o general.o cdata.o
#
sph.o: sph.f90 sph_tree.o sph_kernel.o general.o cdata.o normal_fluid.o output.o forcing.o
#
sph_tree.o: sph_tree.f90 sph_kernel.o general.o cdata.o 
#
sph_kernel.o: sph_kernel.f90 general.o cdata.o 
#
sph_interface.o: sph_interface.f90 sph.o general.o cdata.o normal_fluid.o output.o forcing.o matrix.o
#
hamiltonian.o: hamiltonian.f90 tree.o timestep.o output.o general.o cdata.o
#
