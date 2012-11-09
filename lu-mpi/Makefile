MPI_CC = mpicxx
MPI_CCFLAGS = -DMPICH_IGNORE_CXX_SEEK 

BIN = lu-mpi

all: $(BIN)

lu-mpi: lu-mpi.cpp
	$(MPI_CC) $(MPI_CCFLAGS) -o $@ $<

clean:
	$(RM) $(BIN)

