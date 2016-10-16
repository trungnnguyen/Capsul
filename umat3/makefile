OBJ_DIR=$(HOME)/workspace/abaqus/usrlib
FC=gfortran
FFLAGS= -O3 -fPIC -std=f2003 -ffree-form -fopt-info -fmax-errors=1 -fbounds-check -Werror
SRC=utils.f90	algebra.f90	deformation.f90 hardening.f90 crystal.f90 init.f90 expansion.f90 udisp.f90 umat_cp3.f90 
OBJECTS := $(patsubst %.f90, %-std.o, $(SRC)) capsul.o
capsul: $(OBJECTS)
	rm -fr $(OBJ_DIR)/*
	cp $^ $(OBJ_DIR)
	abaqus make library=$@ directory=$(OBJ_DIR)
utils-std.o: utils.f90
	$(FC) $(FFLAGS) -c $< -o $@
algebra-std.o: algebra.f90 utils-std.o
	$(FC) $(FFLAGS) -c $< -o $@
deformation-std.o: deformation.f90 utils-std.o algebra-std.o
	$(FC) $(FFLAGS) -c $< -o $@
hardening-std.o:	hardening.f90 utils-std.o algebra-std.o
	$(FC) $(FFLAGS) -c $< -o $@
init-std.o:	init.f90 expansion-std.o utils-std.o algebra-std.o crystal-std.o
	$(FC) $(FFLAGS) -c $< -o $@
crystal-std.o: crystal.f90 utils-std.o algebra-std.o hardening-std.o
	$(FC) $(FFLAGS) -c $< -o $@
expansion-std.o: expansion.f90 utils-std.o
	$(FC) $(FFLAGS) -c $< -o $@
udisp-std.o:	udisp.f90
	$(FC) $(FFLAGS) -c $< -o $@
umat_cp3-std.o: umat_cp3.f90 utils-std.o algebra-std.o init-std.o expansion-std.o crystal-std.o \
	             	deformation-std.o hardening-std.o
	$(FC) $(FFLAGS) -c $< -o $@
capsul.o: capsul.f utils-std.o
	$(FC) $(FFLAGS) -c $< -o $@
clean:
	rm *.mod *.o
.PHONY: capsul clean



