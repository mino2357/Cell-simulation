#pragma once

#include <list>
#include "cellToVert.h"

using namespace std;

class cellToVert;

class vertToCell{
public:
	list<int> data;

	cellToVert *CV;

	void show();
	int& operator[](int);

	void push(int i);
	void size();
	void setDataBase(cellToVert*, int, int);
};
