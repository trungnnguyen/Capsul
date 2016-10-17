OBJ_DIR=$(HOME)/workspace/abaqus/usrlib
FC=gfortran
FFLAGS= -O3 -fPIC -std=f2003 -ffree-form -fopt-info -fmax-errors=1 -fbounds-check -Werror
SRC=utils.f90	algebra.f90	 elasticity.f90 slipGeom.f90 slipKinetics.f90 slipHardening.f90 thermoExpansion.f90 \
		crysPlasMat.f90 crysPlasMatPoint.f90 init.f90 
OBJECTS := $(patsubst %.f90, %-std.o, $(SRC)) capsul.o
capsul: $(OBJECTS)
	rm -fr $(OBJ_DIR)/*
	cp $^ $(OBJ_DIR)
	abaqus make library=$@ directory=$(OBJ_DIR)
utils-std.o: utils.f90
	$(FC) $(FFLAGS) -c $< -o $@
algebra-std.o: algebra.f90 utils-std.o
	$(FC) $(FFLAGS) -c $< -o $@
elasticity-std.o:	elasticity.f90 utils-std.o
	$(FC) $(FFLAGS) -c $< -o $@
slipGeom-std.o:	slipGeom.f90 utils-std.o algebra-std.o
	$(FC) $(FFLAGS) -c $< -o $@
slipKinetics-std.o:	slipKinetics.f90 utils-std.o algebra-std.o
	$(FC) $(FFLAGS) -c $< -o $@
slipHardening-std.o:	slipHardening.f90 utils-std.o algebra-std.o
	$(FC) $(FFLAGS) -c $< -o $@
thermoExpansion-std.o:	thermoExpansion.f90 utils-std.o algebra-std.o
	$(FC) $(FFLAGS) -c $< -o $@
crysPlasMat-std.o:	crysPlasMat.f90 elasticity-std.o slipGeom-std.o slipKinetics-std.o	\
		                slipHardening-std.o thermoExpansion-std.o utils-std.o algebra-std.o
	$(FC) $(FFLAGS) -c $< -o $@
crysPlasMatPoint-std.o:	crysPlasMatPoint.f90 crysPlasMat-std.o
	$(FC) $(FFLAGS) -c $< -o $@
init-std.o:	init.f90 crysPlasMat-std.o
	$(FC) $(FFLAGS) -c $< -o $@
capsul.o: capsul.f utils-std.o algebra-std.o init-std.o
	$(FC) $(FFLAGS) -c $< -o $@
clean:
	rm *.mod *.o
.PHONY: capsul clean


