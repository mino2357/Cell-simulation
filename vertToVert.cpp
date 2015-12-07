#include "vertToVert.h"
#include <iostream>
#include <list>
#include <iterator>
#include <cassert>
#include "cellToVert.h"

using namespace std;

vertToVert::vertToVert(){

}

vertToVert::~vertToVert(){

}

//bool isUnique(list<int> data, int n);
// n is my number

void vertToVert::setDataBase(cellToVert *CV, vertToCell *VC, int n){
	buf.clear();
	for(int i=0; i<(int)VC[n].data.size(); i++){
		for(int j=0; j<(int)CV[VC[n][i]].data.size(); j++){
			if(n == CV[VC[n][i]][j]){
				buf.insert(CV[VC[n][i]][j-1]);
				buf.insert(CV[VC[n][i]][j+1]);	
			}
		}
	}
}

void vertToVert::size(){
	buf.size();
}

//void vertToVert::setDataBase(cellToVert *CVdata, int Ncell, int n){
	/*for(int i=0; i<Ncell; i++){
		for(int j=0; j<(int)CVdata[i].data.size(); j++){
			if(n == CVdata[i][j]){
				data.push_back(CVdata[i][j-1]);
				data.push_back(CVdata[i][j+1]);
			}
		}
	}*/
	/*
	list<int>::iterator it = data.begin();

	if(n==0) data.push_back(13);
	if(n==0) data.push_back(4);

	if(n==1) data.push_back(4);
	if(n==1) data.push_back(2);
	
	if(n==2) data.push_back(1);
	if(n==2) data.push_back(3);
	if(n==2) data.push_back(5);
	
	if(n==3) data.push_back(2);
	if(n==3) data.push_back(4);
	if(n==3) data.push_back(11);
	
	if(n==4) data.push_back(3);
	if(n==4) data.push_back(1);
	if(n==4) data.push_back(0);
	
	if(n==5) data.push_back(2);
	if(n==5) data.push_back(6);
	
	if(n==6) data.push_back(5);
	if(n==6) data.push_back(7);
	
	if(n==7) data.push_back(6);
	if(n==7) data.push_back(8);
	
	if(n==8) data.push_back(7);
	if(n==8) data.push_back(11);
	if(n==8) data.push_back(12);
	
	if(n==9) data.push_back(8);
	if(n==9) data.push_back(10);
	
	if(n==10) data.push_back(9);
	if(n==10) data.push_back(11);
	if(n==10) data.push_back(12);
	
	if(n==11) data.push_back(8);
	if(n==11) data.push_back(3);
	if(n==11) data.push_back(10);
	
	if(n==12) data.push_back(10);
	if(n==12) data.push_back(13);
	
	if(n==13) data.push_back(12);
	if(n==13) data.push_back(0);
	*/
	/*it = begin();
	while(it != end()){
		cout << *it << " ";
		it++;
	}
	cout << "\n";
	*/
//}
/*
bool isUnique(list<int> data, int n){
	bool flag = true;
	list<int>::iterator it = data.begin();
	while(it != data.end()){
		if(*it == n){
			return false;
		}
		it++;
	}
	return flag;
}
*/

// for list

/*
void vertToVert::show(){
	list<int>::iterator it = begin();
	while(it != end()){
		cout << *it << " ";
		it++;
	}
	cout << "\n";
}
*/

void vertToVert::show(){
	set<int>::iterator it = buf.begin();
	while(it != buf.end()){
		cout << *it << " ";
		it++;
	}
	cout << "\n";
}

// for list
/*
int& vertToVert::operator[](int i){

	//assert(0 <= i && i < (int)data.size());
	if(i>=0 && i<(int)data.size()){
		list<int>::iterator it = data.begin();
		advance(it, i);
		return *it;
	}

	if(i == -1){
		return data.back();
	}

	if(i == (int)data.size()){
		return data.front();
	}

	cout << "error 2499457" << endl;
	// naosu
	return data.front();
}
*/

int& vertToVert::operator[](int i){

	//assert(0 <= i && i < (int)data.size());
	if(i>=0 && i<(int)buf.size()){
		set<int>::iterator it = buf.begin();
		advance(it, i);
		return (int&)*it;
	}

	if(i == -1){
		return (int&)*buf.end();
	}

	if(i == (int)buf.size()){
		return (int&)*buf.begin();
	}

	cout << "error 2499457" << endl;
	// naosu
	return (int&)*buf.begin();
}


list<int>::iterator vertToVert::begin(){
	return data.begin();
}

list<int>::iterator vertToVert::end(){
	return data.end();
}

