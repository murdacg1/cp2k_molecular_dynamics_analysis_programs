# module load gcc/4.9.2
FC = gfortran

OBJ1 = \
moments_wannier_centers.o

#OP = -g -Wunused -O3 -ffast-math -funroll-loops -ftree-vectorize -march=native -msse3 -ffree-form
OP = -O3 -ffast-math -funroll-loops -ftree-vectorize -march=native -msse3 -ffree-form

BIN1 = moments_wannier_centers.x

# =================================================================

clean:
	rm -f *.o *.mod *.x

# =================================================================

.f.o:
	$(FC) -c $(OP) $<

$(BIN1): $(OBJ1)
	$(FC) $(OBJ1) -o $(BIN1)
