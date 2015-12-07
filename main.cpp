#include <cstdio>
#include <iostream>
#include <list>
#include <iterator>
#include "vertex.h"
#include "cellToVert.h"
#include "vertToVert.h"
#include "update.h"
#include "myglsc.h"
#include "vertToCell.h"
#include "readfile.h"
#include "makePolygon.h"
#include "writeMesh.h"
#include "functions.h"

using namespace std;

const int INTV = 10;

int main(){

	//ファイル名
	string meshFileName = "TT";

	//FreeFem++実行　引数： メッシュファイル名　分割数
	executeFreeFem(meshFileName, 8);

	int Nvert = -1;
	int Ncell = -1;
	
	vertex *vert;
	cellToVert *CV;
	
	//三角形分割データ読み込み　データをvert, CV, Nvert, Ncell に渡す
	readSettingData(meshFileName, vert, CV, &Nvert, &Ncell);
	
	cout << "＝＝＝三角形分割＝＝＝" << endl;
	cout << "頂点数："  << Nvert << " 三角形数" << Ncell << endl;
	
	//メモリ確保
	vertToVert *VV = new vertToVert[Nvert];
	vertToCell *VC = new vertToCell[Nvert];
	
	//データベース設定
	for(int i=0; i<Nvert; i++){
		VC[i].setDataBase(CV, Ncell, i);
	}
	for(int i=0; i<Nvert; i++){
		VV[i].setDataBase(CV, VC, i);
	}
	
	/***************************************************/
	/*                 多角形構成                      */
	/***************************************************/	
	//三角形分割から多角形分割を構成する
	makePolygon mp;
	
	mp.configurePolygon(Nvert, Ncell, VV, vert, CV);

	Ncell = mp.Npolycell;
	Nvert = mp.Npolyvert;
	
	cout << "＝＝＝多角形分割＝＝＝" << endl;
	cout << "頂点数："  << mp.Npolyvert << " 多角形数：" << mp.Npolycell << endl;
	//return 0;
	
	//メモリ確保
	delete[] vert;
	delete[] CV;
	delete[] VV;
	delete[] VC;

	vert = new vertex[Nvert];
	CV = new cellToVert[Ncell];
	VV = new vertToVert[Nvert];
	VC = new vertToCell[Nvert];

	//違反?
	CV = mp.CVpoly;
	vert = mp.v;

	//データベース設定
	for(int i=0; i<Nvert; i++){
		VC[i].setDataBase(CV, Ncell, i);
	}
	for(int i=0; i<Nvert; i++){
		VV[i].setDataBase(CV, VC, i);
	}

	
	/*****************************************************/
	/*                構成終わり                         */
	/*****************************************************/

	//頂点情報のメモリ確保　初期化
	update Update(Nvert, Ncell);

	Update.input(vert, VV, CV, VC);
	Update.setV_std(Ncell);

	//G.L.S.C.3D設定
	myglsc g;
	//g.captureSet();

	bool flag = false;
	bool change = false;

	//初期状態表示
	while(!flag){
		g.clearScreen();
		/* Scele 0 */
		Update.getPotential(Ncell, vert, CV);
		g.sellScale(0);
		g.drawBox(vert, Nvert);
		g.drawXY();
		g.drawPolygon(vert, CV, Ncell);
		
		/* Scele 1 */
		g.sellScale(1);
		g.drawBox2(Update.U);
		g.UI();

		/* Scele 2 */
		g.sellScale(2);
		g.drawBox3(&flag, &change);
		
		//g.rendering();
		g.renderingWithSleep(0.1);
	}

	//cout << __LINE__ << endl;

	flag = true;

	int time = 0;
	//タイムループ
	for(int i=0; ; i++){
		if(flag){
			//Update.Vertex(Nvert, vert, CV, VV, VC);
			Update.Vertex(7, Nvert, vert, CV, VV, VC);
			Update.interchange();
		}

		if(i%INTV == 0){
			Update.getPotential(Ncell, vert, CV);
			
			//データベース設定
			for(int i=0; i<Nvert; i++){
				VC[i].setDataBase(CV, Ncell, i);
			}
			for(int i=0; i<Nvert; i++){
				VV[i].setDataBase(CV, VC, i);
			}
			
			
			//printf("%15.15lf %15.15lf %15.15lf %15.15lf %15.15lf\n", Update.t, Update.dt, Update.U, Update.sigma*Update.S, Update.kappa*Update.V);
			//printf("%15.15lf %15.15lf %15.15lf\n", Update.U, Update.sigma*Update.S, Update.kappa*Update.V);
			g.clearScreen();
			/* Scele 0 */
			g.sellScale(0);
			g.drawBox(vert, Nvert, &flag, &change);
			g.drawXY();
			g.drawPolygon(vert, CV, Ncell);
			
			/* Scele 1 */
			g.sellScale(1);
			g.drawBox2(Update.U);
			g.UI();
			
			/* Scele 2 */
			g.sellScale(2);
			g.drawBox3(&flag, &change);
			
			
			g.renderingWithSleep(0.0);
			//g.rendering();
			//g.capture();
			
			time++;
		}
	}

	delete[] vert;
	delete[] CV;
	delete[] VV;
	delete[] VC;

	return 0;
}
