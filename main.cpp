/*
 * filename: main.cpp
 * author: hadrihilmi@gmail.com
 */

// includes
#include <iostream>

unsigned int N 100 // default 
unsigned int opt 0 // default 0 (sequential)

using namespace std;

// matrix initializationn
void init(int size) {
	float *x = new float[size];
	float *y = new float[size];
}

// main program
int main(int argc, char *argv[]) {

	if(argc < 3) {
		cout << "Usage: $./lu <n> <opt>\n" <<
			"n\tSize of linear system (real number)\n" <<
			"opt\tOperation. 0 - sequential. 1 - openMP. 2 - MPI\n" <<
			"\nExample: $./lu 900 0\n" << endl;
	}

	// catch size of matrix from passing argument
	N = argv[1];

	// catch operation chosen
	opt = argv[2];

	// select engine
	switch(opt) {
		case 0:
			// call sequential engine
			break;
		case 1:
			// call openMP engine
			break;
		case 2:
			// call MPI engine
			break;
	};

	return 0;
}
