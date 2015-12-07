#include <iostream>
#include <iterator>
#include <list>
#include "cellToVert.h"
#include "vertToVert.h"
#include "vertex.h"

using namespace std;

cellToVert::cellToVert(){

}

cellToVert::~cellToVert(){

}

/*
set<int>::iterator cellToVert::begin(){
	return data.begin();
}

set<int>::iterator cellToVert::end(){
	return data.end();
}
*/

void cellToVert::show(){
	list<int>::iterator it = data.begin();
	while(it != data.end()){
		cout << *it << " ";
		it++;
	}
	cout <<"\n";
}

void cellToVert::size(){
	data.size();
}

int& cellToVert::operator[](int i){

	if(i>=0 && i<(int)data.size()){
		list<int>::iterator it = data.begin();
		advance(it, i);
		return *it;
	}

	if(i == -1){
		return data.back();
	}

	if(i == (int)data.size()){
		return data.front();
	}

	cout << "error 2499457" << endl;
	// naosu
	return data.front();
}

double cellToVert::calcVolume(vertex *vert){
	 vertex volume;

	for(int i=0; i<(int)data.size(); i++){
		volume = volume + vert[operator[](i)] % vert[operator[](i+1)];
	}

	return volume.norm()/2.;
}

double cellToVert::calcSurfaceArea(vertex *vert){
	double surfaceArea = 0.;

	for(int i=0; i<(int)data.size(); i++){
		surfaceArea += (vert[operator[](i)] - vert[operator[](i+1)]).norm();
	}

	return surfaceArea;
}

void cellToVert::push(int i){
	data.push_back(i);
}
