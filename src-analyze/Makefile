# FC = gfortran-mp-4.3
FC = gfortran

OBJ1 = \
common_data_mod.o \
binning_energy_mod.o \
reading_mod.o \
minus_com_mod.o \
writing_mod.o \
writing_all_mod.o \
pos_mod.o \
msd_mod.o \
d_msd_mod.o \
vacf_mod.o \
d_vacf_mod.o \
binning_mod.o \
info_mod.o \
charge_mod.o \
dipole_mod.o \
nddo_restraints_mod.o \
analyze.o

OBJ2 = \
d_from_gks_vacf.o

OBJ3 = \
mu.o

OBJ4 = \
mu3.o

OP = -O3 -ffast-math -funroll-loops -ftree-vectorize -march=native -msse3 -ffree-form

BIN1 = analyze.x

BIN2 = d_from_gks_vacf.x

BIN3 = mu.x

BIN4 = mu3.x

# =================================================================

clean:
	rm -f *.o *.mod *.x

# =================================================================

.f.o:
	$(FC) -c $(OP) $<

$(BIN1): $(OBJ1)
	$(FC) $(OBJ1) -o $(BIN1)

$(BIN2): $(OBJ2)
	$(FC) $(OBJ2) -o $(BIN2)

$(BIN3): $(OBJ3)
	$(FC) $(OBJ3) -o $(BIN3)

$(BIN4): $(OBJ4)
	$(FC) $(OBJ4) -o $(BIN4)
