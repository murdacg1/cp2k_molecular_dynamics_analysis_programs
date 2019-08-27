FC = pgf95

OBJ = \
common_data_mod.o \
reading_mod.o \
get_hbond_mod.o \
hbond.o

OP = -O2 -Mfreeform

BIN = hbond.x

# =================================================================

clean:
	rm -f *.o *.mod *.x

# =================================================================

.f.o:
	$(FC) -c $(OP) $<

$(BIN): $(OBJ)
	$(FC) $(OBJ) -o $(BIN)
