OBJ_DIR=$(HOME)/workspace/abaqus/usrlib
FC=gfortran
FFLAGS=-c -O3 -fPIC -std=f2003 -ffree-form -fopt-info -fmax-errors=1 -fbounds-check -Werror
SRC=utils.f90	algebra.f90	deformation.f90 hardening.f90 crystal.f90 init.f90 expansion.f90 umat_cp2.f90 
OBJECTS := $(patsubst %.f90, %-std.o, $(SRC)) capsul.o
capsul: $(OBJECTS)
	rm -fr $(OBJ_DIR)/*
	cp $^ $(OBJ_DIR)
	abaqus make library=$@ directory=$(OBJ_DIR)
utils-std.o: utils.f90
	$(FC) $(FFLAGS) $< -o $@
algebra-std.o: algebra.f90 utils.f90
	$(FC) $(FFLAGS) $< -o $@
deformation-std.o: deformation.f90 utils.f90 algebra.f90
	$(FC) $(FFLAGS) $< -o $@
hardening-std.o:	hardening.f90 utils.f90 algebra.f90
	$(FC) $(FFLAGS) $< -o $@
init-std.o:	init.f90 utils.f90 algebra.f90 crystal.f90 
	$(FC) $(FFLAGS) $< -o $@
crystal-std.o: crystal.f90 utils.f90 algebra.f90 hardening.f90
	$(FC) $(FFLAGS) $< -o $@
expansion-std.o: expansion.f90 utils.f90
	$(FC) $(FFLAGS) $< -o $@
umat_cp2-std.o: umat_cp2.f90 utils.f90 algebra.f90 init.f90 expansion.f90 crystal.f90 deformation.f90 hardening.f90
	$(FC) $(FFLAGS) $< -o $@
capsul.o: capsul.f utils.f90 
	$(FC) $(FFLAGS) $< -o $@
clean:
	rm *.mod *.o
.PHONY: capsul clean



