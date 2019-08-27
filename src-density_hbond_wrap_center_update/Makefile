# module load gcc
FC = gfortran

OBJ = \
common_data_mod.o \
reading_mod.o \
get_density_mod.o \
get_hbond_mod.o \
density_hbond.o

OP = -O3 -ffast-math -funroll-loops -ftree-vectorize -march=native -msse3 -ffree-form

BIN = density_hbond.x

# =================================================================

clean:
	rm -f *.o *.mod *.x

# =================================================================

.f.o:
	$(FC) -c $(OP) $<

$(BIN): $(OBJ)
	$(FC) $(OBJ) -o $(BIN)
