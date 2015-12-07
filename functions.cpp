#include "functions.h"
#include "writeMesh.h"
#include <iostream>
#include <string>
#include <cstdio>
#include <cstdlib>

using namespace std;

void executeFreeFem(string fn, int cut){

	string filename = fn;

	writeMesh wm(filename);
	wm.openWrite(cut);

	string ls = "FreeFem++ " + filename + ".pde";

	system(ls.c_str());

}
