#include "writeMesh.h"
#include <fstream>

using namespace std;

writeMesh::writeMesh(string str){
	filename = str;
}

void writeMesh::openWrite(int cut){
	string name = filename + ".pde";
	ofstream ofs(name.c_str());

	ofs << "border a(t=0, 3.0){x=-1.5+t; y=-1.5;};" << endl;
	ofs << "border b(t=0, 3.0){x=1.5; y=-1.5+t;};" << endl;
	ofs << "border c(t=0, 3.0){x=1.5-t; y=1.5;};" << endl;
	ofs << "border d(t=0, 3.0){x=-1.5; y=1.5-t;};" << endl;
	
	ofs << "int cut = " << cut << ";" << endl;
	ofs << "mesh Th = buildmesh(a(cut)+b(cut)+c(cut)+d(cut));" << endl;

	//ofs << "plot(Th);" << endl;

	ofs << "savemesh(Th, \"" << filename << ".msh\");" << endl;
}
