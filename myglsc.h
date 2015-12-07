#pragma once
#include "vertex.h"
#include "cellToVert.h"

class myglsc{
public:
	//Window Size
	double R;
	int WINDOW_SIZE_X;
	int WINDOW_SIZE_Y;
	
	//def scale
	double x_left;
	double x_right;
	double y_bottom;
	double y_top;

	double x_left_std;
	double y_top_std;

	double width_std;
	double height_std;

	//mouse position
	int downX;
	int downY;
	int x;
	int y;
	
	int changeX;
	int changeY;

	myglsc();
	
	void sellScale(int);
	void clearScreen();
	void testWindow();
	void drawXY();
	
	void drawBox(vertex*, int, bool*, bool*);
	void drawBox(vertex*, int);
	
	void rendering();
	void renderingWithSleep(double);
	void drawPolygon(vertex *, cellToVert *, int);
	void captureSet();
	void capture();

	/* UI */

	void UI();
	
	void drawBox2(double);
	void drawBox3(bool *, bool *);

	void move(vertex*, int);

	double funcX(double);
	double funcY(double);
};
/*
void affiane(vertex *){

}*/
