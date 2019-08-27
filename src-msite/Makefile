# FC = gfortran-mp-4.3
FC = gfortran

OBJ = \
common_data_mod.o \
reading_mod.o \
get_msite_mod.o \
writing_mod.o \
msite.o

# OP = -O3 -ffast-math -funroll-loops -ftree-vectorize -march=native -msse3 -ffree-form
OP = -ffree-form

BIN = msite.x

# =================================================================

clean:
	rm -f *.o *.mod *.x

# =================================================================

.f.o:
	$(FC) -c $(OP) $<

$(BIN): $(OBJ)
	$(FC) $(OBJ) -o $(BIN)
