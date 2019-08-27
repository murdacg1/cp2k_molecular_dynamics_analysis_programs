# FC = gfortran-mp-4.3
FC = gfortran

OBJ = \
combine_histograms.o

OP = -O3 -ffast-math -funroll-loops -ftree-vectorize -msse3 -ffree-form

BIN = combine_histograms.x

# =================================================================

clean:
	rm -f *.o *.mod *.x

# =================================================================

.f.o:
	$(FC) -c $(OP) $<

$(BIN): $(OBJ)
	$(FC) $(OBJ) -o $(BIN)
