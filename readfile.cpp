#include <iostream>
#include <cstdlib>
#include <string>
#include <fstream>
#include "readfile.h"
#include "vertex.h"
#include "cellToVert.h"

using namespace std;

void readSettingData(string filename, vertex *&v, cellToVert *&CV, int *V_N, int *CELL_NUM){

	ifstream fin;

	string str;
	
	string meshfile = filename + ".msh";

	//fin.open("mesh002.msh");
	fin.open(meshfile.c_str());
	if(!fin){
		cout << "error" <<endl;
		exit(1);
	}

	const int maxPoly = 11; //saidai kakkei
	const int DATA_MAX = 50000;

	double a[DATA_MAX][maxPoly];

	int j=0;
	while(getline(fin, str)){
		for(int i=0; i<maxPoly; i++){
			a[j][i] = 0.0;
		}
		sscanf(str.data(), "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &a[j][0], &a[j][1] ,&a[j][2], &a[j][3], &a[j][4], &a[j][5], &a[j][6], &a[j][7], &a[j][8], &a[j][9]);
		j++;
		if(j>=DATA_MAX){
			cout << "error001" <<endl;
			exit(1);
		}
	}

	
	*V_N = a[0][0];
	*CELL_NUM = (int)a[0][1];
	
	int vertex_start = 1;
	//int vertex_end = a[0][0];
	int plgn_start = a[0][0] + 1;
	int plgn_end = plgn_start + a[0][1] - 1;

	v = new vertex[*V_N];
	
	//頂点の設定
	for(int i=0; i<*V_N; i++){
		//v[i].x = 0.2*(a[vertex_start + i][0] - 5);
		//v[i].y = 0.2*(a[vertex_start + i][1] - 5);
		v[i].x = a[vertex_start + i][0];
		v[i].y = a[vertex_start + i][1];
	}

	CV = new cellToVert[*CELL_NUM];

	// CVデータの設定
	int k = 0;
	for(int i=plgn_start; i<=plgn_end; i++){
		//for(int j=1; j<=a[i][0]; j++){
		for(int j=0; j<3; j++){
			CV[k].push(a[i][j]-1);
		}
		k++;
	}
}
