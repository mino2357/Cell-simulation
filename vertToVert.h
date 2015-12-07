#pragma once
#include <list>
#include <set>
#include "cellToVert.h"
#include "vertToCell.h"

class cellToVert;
class vertToCell;

using namespace std;

class vertToVert{
public:
	list<int> data;
	set<int> buf;

	cellToVert *pCVdata;
	vertToVert *pVCdata;

	vertToVert();
	~vertToVert();

	//void setDataBase(cellToVert *,int, int);
	void setDataBase(cellToVert *, vertToCell*, int);

	list<int>::iterator begin();
	list<int>::iterator end();

	void show();
	void size();
	int& operator[](int);
};
