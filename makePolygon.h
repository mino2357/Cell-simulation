#pragma once
#include <vector>
#include "vertToVert.h"
#include "vertex.h"
#include "cellToVert.h"

class vertToVert;
class vertex;
class cellToVert;

class makePolygon{
public:
	int Npolycell;
	int Npolyvert;

	int temp;
	vertex *v;
	cellToVert *CVpoly;

	makePolygon();
	void initCVpoly();
	void findTriangleNumber(int, list<int>::iterator , int, cellToVert *, int);
	int findCV(int , int, int, cellToVert *, int);
	int search(int mv, int new_, int old_, vertToVert *VV);
	void configurePolygon(int Nvert, int Ncell, vertToVert *, vertex *vert, cellToVert *CV);
	void remakePolygon(vertex *, cellToVert *, int *Nvert, int *Ncell);
};
