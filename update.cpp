#include "update.h"
#include "vertex.h"
#include "cellToVert.h"
#include "vertToVert.h"
#include "vertToCell.h"
#include "Embedded_Explicit_Runge_Kutta_Multiple.c"
#include <vector>
#include <set>

using namespace std;

update *update::only_one_pointer = 0;

update::update(int Nvert, int NC){
	
	if(only_one_pointer == 0){

		dS = new vertex[Nvert];
		dV = new vertex[Nvert];

		N_V = Nvert;
		Ncell = NC;

		t = 0.;
		dt = 0.001;
		V_std = 1.1;
		sigma = 0.001;
		kappa = 1.0;

		Rtol = 1.0e-15;
		Atol = 1.0e-15;

		d_min = 1.0e-3;

		S = 0.;
		V = 0.;
		U = 0.;

		only_one_pointer = this;
	
	}else{
		return;
	}
}

update::~update(){
	delete[] dS;
	delete[] dV;
}

void update::interchange(){
	getThresholdSet();

	list<pair<int, int> >::iterator itv = thresholdVertex.begin();
	list<pair<int, int> >::iterator itc = thresholdCell.begin();
	list<pair<int, int> >::iterator ita = anotherCell.begin();

	while(itv != thresholdVertex.end()){
		//cout << (*itv).first << " " << (*itv).second << endl;
		//cout << (*itc).first << " " << (*itc).second << endl;
		//cout << (*ita).first << " " << (*ita).second << endl;
		
		/*     頂点の繋ぎ替え     */

		// remove
		CV[(*itc).first].data.remove((*itv).first);
		CV[(*itc).second].data.remove((*itv).second);

		//add
		
		sort((*ita).first, (*itv).second);
		sort((*ita).second, (*itv).first);

		CV[(*ita).first].data.push_front((*itv).first);
		CV[(*ita).second].data.push_front((*itv).second);

		//CV[(*ita).first].data.insert()

		/**************************/
		
		itv++;
		itc++;
		ita++;
	}
}

void update::sort(int cn, int vn){
	//cout << cn << " " << vn << endl;

	vector<int> data;

	int dataSize = (int)CV[cn].data.size();

	for(int i=0; i<dataSize; i++){
		data.push_back(CV[cn][i]);
	}

	//vn ha nannbannme ni aruka?
	int num = -1;

	for(int i=0; i<dataSize; i++){
		if(vn == data[i]){
			num = i;
		}
	}

	// clear!
	CV[cn].data.clear();

	//vn senntou narabikae
	for(int i=num; i<dataSize; i++){
		CV[cn].data.push_back(data[i]);
	}

	for(int i=0; i<num; i++){
		CV[cn].data.push_back(data[i]);
	}

/*
	for(int i=0; i<CV[cn].data.size(); i++){
		cout << CV[cn][i];
	}
	cout << "\n";
*/
}

void update::getThresholdSet(){
	//初期化
	thresholdVertex.clear();
	thresholdCell.clear();
	anotherCell.clear();

	for(int i=0; i<N_V; i++){
		for(int j=0; j<=i; j++){
	//cout << __LINE__ << endl;
			if(VVlength(i,j) < d_min){
				thresholdVertex.push_back(make_pair(j,i));
				getCell(i, j);
			}
		}
	}
}

void update::getCell(int m, int n){
	int c0 = -1;
	int c1 = -1;

	int time = 0;

	for(int i=0; i<Ncell; i++){
		for(int j=0; j<(int)CV[i].data.size(); j++){
			if((CV[i][j] == m && CV[i][j+1] == n) || (CV[i][j] == n && CV[i][j+1] == m)){
				if(time > 0){
					c1 = i;
					break;
				}
				c0 = i;
				time++;
			}
		}
	}

	vector<int> data;

	//cout << m << " " << n << endl;
	//cout << c0 << " " << c1 << endl;

	int i = m;
	for(int j=0; j<(int)VC[i].data.size(); j++){
		if(VC[i][j] != c0 && VC[i][j] != c1){
			data.push_back(VC[i][j]);
			//cout << VC[i][j] << endl;
		}
	}
	
	i = n;
	for(int j=0; j<(int)VC[i].data.size(); j++){
		if(VC[i][j] != c0 && VC[i][j] != c1){
			data.push_back(VC[i][j]);
			//cout << VC[i][j] << endl;
		}
	}

	//cout << data[0] << " " << data[1] << endl;
	
	thresholdCell.push_back(make_pair(c0, c1));
	anotherCell.push_back(make_pair(data[0], data[1]));
}

double update::VVlength(int n, int m){
	for(int i=0; i<(int)VV[n].buf.size(); i++){
		if(VV[n][i] == m){
			return (v[n] - v[VV[n][i]]).norm();
		}
	}
	return 1000000;
}

void update::setV_std(int Ncell){
	double V=0.;
	for(int i=0; i<Ncell; i++){
		V += CV[i].calcVolume(v);
	}
	V_std = V/Ncell;
}

void update::input(vertex *vert, vertToVert *VVdata, cellToVert *CVdata, vertToCell *VCdata){
	v = vert;
	VV = VVdata;
	CV = CVdata;
	VC = VCdata;
}

//埋め込み型ルンゲ・クッタ法
void update::Vertex(int solverNum, int N_V, vertex *vert, cellToVert *CVdata, vertToVert *VVdata, vertToCell *VCdata){
	input(vert, VVdata, CVdata, VCdata);

	double *u = new double[2 * N_V];

	int num = 0;

	for(int i = 0; i < 2*N_V; i=i+2){
		u[i] = v[num].x;
		u[i+1] = v[num].y;
		num++;
	}


	Multiple_n_dim(solverNum, 2*N_V, u, &t, &dt, Atol, Rtol, func);
	
	/*if(dt>0.1){
		dt = 0.01;
	}*/

	num = 0;
	for(int i=0; i<2*N_V; i=i+2){
		v[num].x = u[i];
		v[num].y = u[i+1];
		num++;
	}

	delete[] u;
}


void update::func(int n, double *u, double *u_dot, double time){

	int count = 0;
	
	update *p = only_one_pointer;
	int num = 0;
	for(int i = 0; i < 2*p->N_V; i=i+2){
		p->v[num].x = u[i];
		p->v[num].y = u[i+1];
		num++;
	}

	p->calcVariationsS(p->N_V, p->v, p->VV);
	p->calcVariationsPlusV(p->N_V, p->v, p->CV, p->VC);

	for(int i=0; i<p->N_V; i++){
		u_dot[count] = - (p->sigma * p->dS[i] + p->kappa * p->dV[i]).x;
		count++;
		u_dot[count] = - (p->sigma * p->dS[i] + p->kappa * p->dV[i]).y;
		count++;
	}
	num = 0;
	for(int i=0; i<2*p->N_V; i=i+2){
		p->v[num].x = u[i];
		p->v[num].y = u[i+1];
		num++;
	}
}

//オイラー法
void update::Vertex(int Nvert, vertex *vert, cellToVert *CVdata, vertToVert *VVdata, vertToCell *VCdata){
	calcVariationsS(Nvert, vert, VVdata);
	calcVariationsPlusV(Nvert, vert, CVdata, VCdata);

	for(int i=0; i<Nvert; i++){
		vert[i] = vert[i] - dt * ( sigma * dS[i] + kappa * dV[i] );
	}
	t += dt;
}

void update::calcVariationsS(int Nvert, vertex *vert, vertToVert *VVdata){
	vertex init;

	for(int i=0; i<Nvert; i++){
		dS[i] = init;
	}

	//vertex v;
	for(int i=0; i<Nvert; i++){
		//list<int>::iterator it = VVdata[i].begin();
		for(int j=0; j<(int)VVdata[i].buf.size(); j++){
			//dS[i] = dS[i] + (vert[i] - vert[*it])/(vert[i] - vert[*it]).norm();
			dS[i] = dS[i] + (vert[i] - vert[VVdata[i][j]])/(vert[i] - vert[VVdata[i][j]]).norm();
			//it++;
		}
	}
}


void update::calcVariationsPlusV(int Nvert, vertex *vert, cellToVert *CV, vertToCell* VC){
	//となりのふたつ
	int a = 0;
	int b = 0;

	vertex init;
	vertex v = init;
	double volume = 0.;

	for(int vv = 0; vv < Nvert; vv++){
		dV[vv] = init;
		//v頂点のまわりの細胞を調査
		for(int i=0; i<(int)VC[vv].data.size(); i++){
			v = init;
			volume = 0.;
			//i番細胞の中を回る
			for(int j=0; j<(int)CV[VC[vv][i]].data.size(); j++){
				v = v + (vert[CV[VC[vv][i]][j]] % vert[CV[VC[vv][i]][j+1]]);
				if(vv == CV[VC[vv][i]][j]){
					a = CV[VC[vv][i]][j-1];
					b = CV[VC[vv][i]][j+1];
				}
			}
			volume = CV[VC[vv][i]].calcVolume(vert);
			dV[vv] = dV[vv] + (volume - V_std) * (v % (vert[a] - vert[b]))/v.norm();
		}
	}
}

void update::getPotential(int Ncell, vertex *v, cellToVert *CV){
	S = 0.;
	V = 0.;
	U = 0.;

	for(int i=0; i<Ncell; i++){
		S += CV[i].calcSurfaceArea(v);
		V += (CV[i].calcVolume(v) - V_std) * (CV[i].calcVolume(v) - V_std);
	}

	U = sigma * S + kappa * V;
}

	// vv vertex  vv点のとなりを表示
	/*
	int vv = 3;
	for(int i=0; i<(int)VC[vv].data.size(); i++){
		for(int j=0; j<(int)CV[VC[vv][i]].data.size(); j++){
			if(vv == CV[VC[vv][i]][j]){
				cout << CV[VC[vv][i]][j-1] << endl;	
				cout << CV[VC[vv][i]][j+1] << endl;	
			}
		}
		cout << "\n";
	}
	cout << "\n";
	*/
	
	/*
	//となりのふたつ
	int a = 0;
	int b = 0;
	
	double volume = 0.;

	int vv = 3;
	dV[vv] = init;
	//v頂点のまわりの細胞を調査
	for(int i=0; i<(int)VC[vv].data.size(); i++){
		v = init;
		volume = 0.;
		//i番細胞の中を回る
		for(int j=0; j<(int)CV[VC[vv][i]].data.size(); j++){
			v = v + (vert[CV[VC[vv][i]][j]] % vert[CV[VC[vv][i]][j+1]]);
			if(vv == CV[VC[vv][i]][j]){
				a = CV[VC[vv][i]][j-1];
				b = CV[VC[vv][i]][j+1];
			}
		}
		volume = CV[VC[vv][i]].calcVolume(vert);
		dV[vv] = dV[vv] + 0.5 * (volume - V_std) * (v % (vert[b] - vert[a]))/v.norm();
	}
	cout << "\n";
	
	vertex init;
	vertex v;
	
	for(int i=0; i<Nvert; i++){
		dV[i] = init;
	}
	*/

