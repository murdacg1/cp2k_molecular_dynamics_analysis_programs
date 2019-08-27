# module load gcc/4.8.2
FC = gfortran

OBJ1 = \
common_data_mod.o \
reading_writing_mod.o \
center_and_wrap.o

#OP = -g -Wunused -O3 -ffast-math -funroll-loops -ftree-vectorize -march=native -msse3 -ffree-form
OP = -O3 -ffast-math -funroll-loops -ftree-vectorize -march=native -msse3 -ffree-form

BIN1 = center_and_wrap.x

# =================================================================

clean:
	rm -f *.o *.mod *.x

# =================================================================

.f.o:
	$(FC) -c $(OP) $<

$(BIN1): $(OBJ1)
	$(FC) $(OBJ1) -o $(BIN1)
