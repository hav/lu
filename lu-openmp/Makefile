CC = g++
CC_FLAGS = -fopenmp

BIN = lu-omp

all: $(BIN)

lu-omp: lu-omp.cpp
	$(CC) $(CC_FLAGS) -o $@ $<

clean:
	$(RM) $(BIN)