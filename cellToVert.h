#pragma once
#include <list>
#include "vertex.h"
#include "vertToVert.h"

class vertToVert;

using namespace std;

class cellToVert{
public:
	list<int> data;

	vertToVert *pVVdata;

	cellToVert();
	~cellToVert();
	
	list<int>::iterator begin();
	list<int>::iterator end();

	void show();
	void size();
	int& operator[](int);

	double calcVolume(vertex *);
	double calcSurfaceArea(vertex *);

	void push(int);
};
