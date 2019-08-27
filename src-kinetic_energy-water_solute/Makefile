# module load gcc
FC = gfortran

OBJ1 = \
common_data_mod.o \
reading_and_calc_ke_mod.o \
binning_mod.o \
kinetic_energy.o

#OP = -g -fbounds-check -Warray-bounds -Wunused -O3 -ffast-math -funroll-loops -ftree-vectorize -march=native -msse3 -ffree-form
OP = -O3 -ffast-math -funroll-loops -ftree-vectorize -march=native -msse3 -ffree-form

BIN1 = kinetic_energy.x

# =================================================================

clean:
	rm -f *.o *.mod *.x

# =================================================================

.f.o:
	$(FC) -c $(OP) $<

$(BIN1): $(OBJ1)
	$(FC) $(OBJ1) -o $(BIN1)
