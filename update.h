#pragma once
#include <iostream>
#include "vertToVert.h"
#include "cellToVert.h"
#include "vertex.h"
#include "vertToCell.h"
#include <list>
#include <utility>

class vertex;
class vertToVert;
class cellToVert;
class vertToCell;

using namespace std;

class update{
	static update *only_one_pointer;
public:
	update(int, int);
	~update();

	/*******************/
	vertex *v;
	vertToVert *VV;
	cellToVert *CV;
	vertToCell *VC;

	void input(vertex *, vertToVert *, cellToVert *, vertToCell *);
	/*******************/

	//各定数
	double t;
	double dt;
	double V_std;
	double sigma;
	double kappa;

	double Rtol;
	double Atol;

	double d_min;

	int N_V;
	int Ncell;

	list<pair<int, int> > thresholdVertex;
	list<pair<int, int> > thresholdCell;
	list<pair<int, int> > anotherCell;

	//全体の表面積，体積
	double S;
	double V;
	double U;

	vertex *dS;
	vertex *dV;

	void Vertex(int Nvert, vertex *, cellToVert*, vertToVert*, vertToCell*);
	void Vertex(int solverNum, int N_V, vertex *, cellToVert*, vertToVert*, vertToCell*);
	void calcVariationsS(int Nvert, vertex *, vertToVert *);
	void calcVariationsPlusV(int Nvert, vertex *, cellToVert *, vertToCell*);

	void setV_std(int Ncell);

	static void func(int n, double *u, double *u_dot, double time);

	void getPotential(int Ncell, vertex *, cellToVert*);

	void getThresholdSet();
	void getCell(int, int);

	void sort(int cn, int vn);

	double VVlength(int, int);
	void interchange();
};
