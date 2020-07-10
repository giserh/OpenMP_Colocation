#include <iostream>
#include <stdlib.h>
//#include <cuda_runtime_api.h>
#include "colocationFinder.h"
using namespace std;


int main(int argc, char**argv) {
	colocationFinder* oColocationFinder = NULL;
	oColocationFinder = new colocationFinder();
	oColocationFinder->Begin(argc,argv);
	cout << "End";

	return 0;
}