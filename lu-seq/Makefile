CC = g++
CC_FLAGS = -Wall -O2 -g

BIN = lu-seq

all: $(BIN)

lu-seq: lu-seq.cpp
	$(CC) $(CC_FLAGS) -o $@ $<

clean:
	$(RM) $(BIN)