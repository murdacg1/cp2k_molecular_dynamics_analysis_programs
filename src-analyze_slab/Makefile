# FC = gfortran-mp-4.3
FC = gfortran

OBJ1 = \
common_data_mod.o \
reading_mod.o \
minus_com_mod.o \
writing_all_mod.o \
analyze.o

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
