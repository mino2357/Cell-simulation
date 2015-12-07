#include <iostream>
#include <string>

using namespace std;

class writeMesh{
public:
	writeMesh(string);
	
	string filename;

	void openWrite(int cut);
};
