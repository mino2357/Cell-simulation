#include "myglsc.h"
#include "glsc3d.h"
#include <string>
#include <sstream>
#include <cmath>

using namespace std;


myglsc::myglsc(){
	//Window Size
	R = 1.0;
	WINDOW_SIZE_X	 = (int)900*R;
	WINDOW_SIZE_Y  = (int)700*R;
	
    g_init_core("Graph",WINDOW_SIZE_X,WINDOW_SIZE_Y,0,0,1,1,1,0,0,0);
	
	//def scale
	x_left = -2.0;
	x_right = 2.0;
	y_bottom = -2.0;
	y_top = 2.0;

	x_left_std = 20.0;
	y_top_std = 20.0;
	
	width_std = 500*R;
	height_std = 500*R;
	
	g_def_scale_2D(
		0,
		x_left, x_right,
		y_bottom, y_top,
		x_left_std, y_top_std,
		width_std, height_std);
	
	//def scale 1
	x_left = -2.0;
	x_right = 2.0;
	y_bottom = -2.0;
	y_top = 2.0;

	x_left_std = 40.0 + width_std;
	y_top_std = 20.0;
	
	width_std = 340*R;
	height_std = 500*R;
	
	g_def_scale_2D(
		1,
		x_left, x_right,
		y_bottom, y_top,
		x_left_std, y_top_std,
		width_std, height_std);
	
	//def scale 2
	x_left_std = 20.0;
	y_top_std = 40.0 + 500*R;
	
	width_std = 860*R;
	height_std = 140*R;
	
	x_left = 20;
	x_right = WINDOW_SIZE_X-20;
	y_bottom = WINDOW_SIZE_Y-20;
	y_top = 20 + 500*R + 20;

/*
	x_left = 20.;
	x_right = WINDOW_SIZE_X;
	y_bottom = 0.;
	y_top = WINDOW_SIZE_Y;
*/
	g_def_scale_2D(
		2,
		x_left, x_right,
		y_bottom, y_top,
		x_left_std, y_top_std,
		width_std, height_std);
	
	//def scale
	x_left = -2.0;
	x_right = 2.0;
	y_bottom = -2.0;
	y_top = 2.0;

	//red
	g_def_line(0, 1.0, 0.0, 0.0, 1.0, 0, 0);
	
	//green
	g_def_line(1, 0.0, 1.0, 0.0, 1.0, 0, 0);
	
	//blue
	g_def_line(2, 0.0, 0.0, 1.0, 1.0, 0, 0);
	
	g_sel_scale_2D(0);
	
	g_def_text(0, 0, 0, 0, 1, 4);
	g_def_text(1, 1, 0, 0, 1, 4);
	g_def_text(2, 0, 1, 0, 1, 4);
	g_def_text(3, 0, 0, 1, 1, 4);
	//g_area_color_2D(0.5, 0.5, 0.5, 0.5);
	
	
	//mouse position
	
	downX = -10000;
	downY = -10000;
	
	x = -100000;
	y = -100000;
	
	changeX = -100000;
	changeY = -100000;
	
	x_left_std = 20.0;
	y_top_std = 20.0;
	
	width_std = 500*R;
	height_std = 500*R;
	
	//def scale
	x_left = -2.0;
	x_right = 2.0;
	y_bottom = -2.0;
	y_top = 2.0;
}

double myglsc::funcX(double X){
	double A = (X - x_left_std) / (width_std);// +  x_left_std);
	return x_left + A * (x_right - x_left);
}

double myglsc::funcY(double Y){
	double B = (Y - y_top_std) / (height_std);// - y_top_std);
	return -(y_bottom + B * (y_top - y_bottom));
}

void myglsc::sellScale(int i){
	g_sel_scale_2D(i);
}

void myglsc::testWindow(){
	g_finish();
	g_sleep(2);
}

void myglsc::drawXY(){
	g_sel_scale_2D(0);

	g_sel_line(0);
	g_move_2D(-1000.0, 0.0);
	g_plot_2D(1000.0, 0.0);

	g_sel_line(1);
	g_move_2D(0.0, -1000.0);
	g_plot_2D(0.0, 1000.0);

/*
	g_text_standard(1,1,"Hello!");

	g_sel_text(1);
	g_text_2D_virtual(1, 1, "Hello");
	g_sel_text(0);
	g_text_2D_virtual(1, 1, "Hello");
*/}

void myglsc::clearScreen(){
	g_cls();
}

void myglsc::rendering(){
	g_finish();
}

void myglsc::renderingWithSleep(double time){
	g_finish();
	g_sleep(time);
}

void myglsc::drawBox(vertex *v, int Nvert){
	g_line_color(0, 0, 0, 1);
	g_line_width(6);
	g_area_color_2D(1, 0, 0, 1);

	g_box_2D((x_right + x_left)/2., (y_top + y_bottom)/2., x_right-x_left, y_top-y_bottom, (G_WIREFILL)0);
	g_line_width(1);

}

void myglsc::move(vertex *v, int Nvert){
	bool in = false;
	
	if(g_input_state(G_MOUSE_LEFT, &x, &y) == G_DOWN){
		if(x<520 && x>20 && y>20 && y<520){
			changeX = x;
			changeY = y;
			in = true;
		}
	}


	if(in){
		vertex mouse;
		mouse.x = funcX(changeX);
		mouse.y = funcY(changeY);
	
		//cout << mouse.x << " " << mouse.y << endl;

		int cNum = -1;

		double length = 100000;

		for(int i=0; i<Nvert; i++){
			if((v[i] - mouse).norm() < length){
				cNum = i;
				length = (v[i] - mouse).norm();
			}
		}

		v[cNum] = mouse;

		//cout << cNum << endl;

		if(cNum == -1){
			return;
		}
	}
}

void myglsc::drawBox(vertex *v, int Nvert, bool *flag, bool *change){
	g_line_color(0, 0, 0, 1);
	g_line_width(6);
	g_area_color_2D(1, 0, 0, 1);

	g_box_2D((x_right + x_left)/2., (y_top + y_bottom)/2., x_right-x_left, y_top-y_bottom, (G_WIREFILL)0);
	g_line_width(1);

	if(!*flag && *change){
		move(v, Nvert);
	}
}

void myglsc::drawPolygon(vertex *vert, cellToVert *CVdata, int Ncell){
	g_sel_scale_2D(0);
	g_sel_line(2);
	g_line_width(2);

	for(int i=0; i<Ncell; i++){
		g_move_2D(vert[CVdata[i][0]].x, vert[CVdata[i][0]].y);
		for(int j=1; j<=(int)CVdata[i].data.size(); j++){
			g_plot_2D(vert[CVdata[i][j]].x, vert[CVdata[i][j]].y);
		}
	}
	g_line_width(1);
}

void myglsc::captureSet(){
	g_capture_set("");
}

void myglsc::capture(){
	g_capture();
}

void myglsc::UI(){
//	g_sel_text(1);
//	g_text_2D_virtual(100, 68, "Exit");
//	g_text_2D_virtual(100, 93, "Program");
}

void myglsc::drawBox2(double U){
	g_line_color(0, 1, 0, 1);
	g_line_width(2);
	g_area_color_2D(1, 0, 0, 1);

	g_box_2D((x_right + x_left)/2., (y_top + y_bottom)/2., x_right-x_left, y_top-y_bottom, (G_WIREFILL)0);
	g_line_width(1);
	
	g_sel_text(0);

	stringstream ss;

	ss << U;
	//ss.precision(2); 
	//ss.fixed; 
	string data = "U: " + ss.str();

	g_text_2D_virtual(0, 0, data.c_str());

}

void myglsc::drawBox3(bool *flag, bool *change){
	g_line_color(0, 0, 1, 1);
	g_line_width(2);
	g_area_color_2D(1, 0, 0, 1);

	g_box_2D(WINDOW_SIZE_X/2., (WINDOW_SIZE_Y + 500*R + 20)/2., WINDOW_SIZE_Y+160, 140*R, (G_WIREFILL)0);
	g_line_width(1);


	if(g_input_state(G_MOUSE_LEFT, &x, &y) == G_DOWN){
		downX = x;
		downY = y;
		//cout << x << " " << y << endl;
		//g_sleep(1.);
	}
	g_circle_2D(downX, downY, 10, (G_WIREFILL)1);

	/* create box */
	g_line_color(0, 0, 1, 1);
	g_line_width(2);
	g_area_color_2D(0, 0, 1, 0.5);
	//g_box_2D(20+170*1, 540+70, 120, 40, (G_WIREFILL)1);
	//g_box_2D(20+170*2, 540+70, 120, 40, (G_WIREFILL)1);
	//g_box_2D(20+170*3, 540+70, 120, 40, (G_WIREFILL)1);
	//g_box_2D(20+170*4, 540+70, 120, 40, (G_WIREFILL)1);
	
	g_line_color(0, 0, 0, 1);
	g_box_2D(20+170*1, 540+70, 120, 40, (G_WIREFILL)0);
	g_box_2D(20+170*2, 540+70, 120, 40, (G_WIREFILL)0);
	g_box_2D(20+170*3, 540+110,120, 40, (G_WIREFILL)0);
	g_box_2D(20+170*3, 540+30, 120, 40, (G_WIREFILL)0);
	g_box_2D(20+170*4, 540+70, 120, 40, (G_WIREFILL)0);
	g_line_width(1);

	g_sel_text(0);

	//start
	g_text_2D_virtual(170*1, 540+75, "Start");
	if(downX>(20+170*1-60) && downX<(20+170*1+60)
		&& downY<(540+70+20) && downY>(540+70-20)){
		//cout << *flag << endl;
		if(!*flag){
			*flag = true;
		}
		//cout << *flag << endl;
		downX = downY = 100000;
	}
	//stop
	g_text_2D_virtual(170*2, 540+75, "Stop");
	if(downX>(20+170*2-60) && downX<(20+170*2+60)
		&& downY<(540+70+20) && downY>(540+70-20)){
		//cout << *flag << endl;
		if(*flag){
			*flag = false;
		}
		//cout << *flag << endl;
		downX = downY = 100000;
	}
	//change
	g_text_2D_virtual(170*3-10, 540+75-40, "Change");
	if(downX>(20+170*3-60) && downX<(20+170*3+60)
		&& downY<(540+70+20-40) && downY>(540+70-20-40)){
		if(!*change){
			*change = true;
		}
		downX = downY = 100000;
	}
	//save
	g_text_2D_virtual(170*3, 540+75+40, "Save");
	if(downX>(20+170*3-60) && downX<(20+170*3+60)
		&& downY<(540+70+20+40) && downY>(540+70-20+40)){
		if(*change){
			*change = false;
		}
		downX = downY = 100000;
	}
	//exit
	g_text_2D_virtual(170*4, 540+75, "Exit");
	if(downX>(20+170*4-60) && downX<(20+170*4+60)
		&& downY<(540+70+20) && downY>(540+70-20)){
		downX = downY = 100000;
		exit(0);	
	}


}
