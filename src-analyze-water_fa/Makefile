# module load gcc
FC = gfortran

OBJ1 = \
common_data_mod.o \
binning_energy_mod.o \
reading_mod.o \
writing_all_mod.o \
pos_mod.o \
pos_sol_mod.o \
binning_mod.o \
dist_mod.o \
writing_dist_mod.o \
analyze.o

#OP = -g -Wunused -O3 -ffast-math -funroll-loops -ftree-vectorize -march=native -msse3 -ffree-form
OP = -O3 -ffast-math -funroll-loops -ftree-vectorize -march=native -msse3 -ffree-form

BIN1 = analyze.x

# =================================================================

clean:
	rm -f *.o *.mod *.x

# =================================================================

.f.o:
	$(FC) -c $(OP) $<

$(BIN1): $(OBJ1)
	$(FC) $(OBJ1) -o $(BIN1)
