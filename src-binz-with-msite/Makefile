# FC = gfortran-mp-4.3
FC = gfortran

OBJ = \
common_data_mod.o \
reading_mod.o \
binning_z_mod.o \
binz.o

OP = -O3 -ffast-math -funroll-loops -ftree-vectorize -msse3 -ffree-form

BIN = binz.x

# =================================================================

clean:
	rm -f *.o *.mod *.x

# =================================================================

.f.o:
	$(FC) -c $(OP) $<

$(BIN): $(OBJ)
	$(FC) $(OBJ) -o $(BIN)
