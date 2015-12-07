#pragma once

#include <iostream>
#include <cmath>
using namespace std;

class vertex{
public:
	double x;
	double y;
	double z;

	int b; //border -> -1

	vertex();
	vertex(double, double);
	vertex(double, double, double);
	void set(double, double);
	vertex operator+(vertex obj);
	vertex operator-(vertex obj);
	vertex operator-();
	vertex operator/(double a);
	double operator*(vertex obj);
	vertex operator%(vertex obj);
	friend vertex operator*(double, vertex obj);
	void show();
	double norm();
	void rot(double);
};
