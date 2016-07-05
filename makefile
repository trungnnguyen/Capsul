OBJ_DIR=$(HOME)/workspace/abaqus/usrlib
FC=gfortran
FFLAGS=-c -O3 -fPIC -std=f2003 -ffree-form -fopt-info -fmax-errors=1 -fbounds-check -Werror
MAIN_SRC=capsul.f
COMMON_SRC=utils.f90	algebra.f90	deformation.f90 hardening.f90 init.f90
OBJECTS := $(patsubst %.f90, %-std.o, $(COMMON_SRC))
UMAT_SRC=./umat
capsul:$(OBJECTS)
	make -C $(COMMON_SRC)
	make -C $(UMAT_SRC)
	abaqus make library=${MAIN%.*} directory=$(DIR)
utils-std.o: utils.f90
	$(FC) $(FFLAGS) $< -o $@
algebra-std.o: algebra.f90 utils.f90
	$(FC) $(FFLAGS) $< -o $@
deformation-std.o: deformation.f90 utils.f90 algebra.f90
	$(FC) $(FFLAGS) $< -o $@
hardening-std.o:	hardening.f90 utils.f90 algebra.f90
	$(FC) $(FFLAGS) $< -o $@
init-std.f90:	init.f90 utils.f90 algebra.f90 crystal.f90 
	$(FC) $(FFLAGS) $< -o $@
clean:
	rm *.mod *.o
.PHONY: capsul clean



