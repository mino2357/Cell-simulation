#include "makePolygon.h"

makePolygon::makePolygon(){
	Npolycell = -1;
	Npolyvert = -1;

	temp = 0;
}

int makePolygon::search(int mv, int new_, int old_, vertToVert *VV){
	for(int i=0; i<(int)VV[mv].buf.size(); i++){
		for(int j=0; j<(int)VV[new_].buf.size(); j++){
			if(VV[mv][i] == VV[new_][j] && old_ != VV[new_][j]){
				return VV[mv][i];
			}
		}
	}
	return -1;
}

void makePolygon::configurePolygon(int Nvert, int Ncell, vertToVert *VV, vertex *vert, cellToVert *CV){
	list<int> gData;
	list<int> centerData;

	/*　構成開始　*/

	int NpolyCount = 0;

	for(int mv = 0; mv<Nvert; mv++){
	//for(int mv = 1; mv<2; mv++){
		int first = VV[mv][0];

		int new_ = first; 
		int old_ = mv;
		int Npoly = 0;
		gData.push_back(first);
		for(int i=0; i<(int)VV[mv].buf.size(); i++){
			int next = search(mv, new_, old_, VV);

			if(next == -1){
				break;
			}else{
				old_ = new_; 
				new_ = next;
				Npoly++;
				if(first == next && Npoly>2){
					//cout << mv << endl;
					//for(list<int>::iterator itr = gData.begin(); itr != gData.end(); ++itr) {
				      //  cout << *itr << " ";	// イテレータの指す先のデータを表示
					//}
					//cout << "\n";
					//cout << findTriangleNumber(mv, gData.begin(), Npoly, CV, Ncell) << endl;
					centerData.push_back(mv);
					gData.clear();
					NpolyCount++;
				}
				gData.push_back(next);
			}
		}
		gData.clear();
	}

	//多角形の個数
	Npolycell = NpolyCount;
	//

	/*　構成終終了　*/

	//cout << "多角形の個数：　" << Npolycell << endl;

	/* 重心構成 */

	v = new vertex[Ncell];
	vertex init;

	for(int i=0; i<Ncell; i++){
		v[i] = init;
	}
	
	for(int i=0; i<Ncell; i++){
		for(int j=0; j<(int)CV[i].data.size(); j++){
			v[i] = v[i] + vert[CV[i][j]];
		}
		v[i] = v[i]/3.;
	}

	Npolyvert = Ncell;

	/* 再構成　代入 */
	gData.clear();
	centerData.clear();

	NpolyCount = 0;
	
	initCVpoly();
	
	for(int mv = 0; mv<Nvert; mv++){
		int first = VV[mv][0];

		int new_ = first; 
		int old_ = mv;
		int Npoly = 0;
		gData.push_back(first);
		for(int i=0; i<(int)VV[mv].buf.size(); i++){
			int next = search(mv, new_, old_, VV);

			if(next == -1){
				break;
			}else{
				old_ = new_; 
				new_ = next;
				Npoly++;
				if(first == next && Npoly>2){
					//cout << mv << endl;
					//for(list<int>::iterator itr = gData.begin(); itr != gData.end(); ++itr) {
				      //  cout << *itr << " ";	// イテレータの指す先のデータを表示
					//}
					//cout << "\n";
					findTriangleNumber(mv, gData.begin(), Npoly, CV, Ncell);
					centerData.push_back(mv);
					gData.clear();
					NpolyCount++;
				}
				gData.push_back(next);
			}
		}
		gData.clear();
	}
}

void makePolygon::findTriangleNumber(int mv, list<int>::iterator it, int size, cellToVert *CV, int Ncell){

	int begin = *it;
	int end = -1;

	//cout << "size : " << size << endl;
	
	for(int i=1; i<size; i++){
		list<int>::iterator ita = it;
		list<int>::iterator itb = ++it;

		//cout << *ita << " "  << *itb << " " << endl;
		//cout << findCV(mv, *ita, *itb, CV, Ncell) << endl;
		end = *itb;

		CVpoly[temp].data.push_back(findCV(mv, *ita, *itb, CV, Ncell));
	}

	//cout << begin << " " << end << endl;

	//cout << findCV(mv, begin, end, CV, Ncell) << endl;;
	CVpoly[temp].data.push_back(findCV(mv, begin, end, CV, Ncell));

	temp++;
}

int makePolygon::findCV(int a, int b, int c, cellToVert *CV, int Ncell){
	for(int i=0; i<Ncell; i++){
		if(      CV[i][0] == a && CV[i][1] == b && CV[i][2] == c){
			return i;
		}else if(CV[i][0] == a && CV[i][2] == b && CV[i][1] == c){
			return i;
		}else if(CV[i][1] == a && CV[i][0] == b && CV[i][2] == c){
			return i;
		}else if(CV[i][1] == a && CV[i][2] == b && CV[i][0] == c){
			return i;
		}else if(CV[i][2] == a && CV[i][0] == b && CV[i][1] == c){
			return i;
		}else if(CV[i][2] == a && CV[i][1] == b && CV[i][0] == c){
			return i;
		}
	}
	cout << "error y8vrneiovoei" << endl;
	return -1;
}

void makePolygon::initCVpoly(){
	CVpoly = new cellToVert[Npolycell];
}
