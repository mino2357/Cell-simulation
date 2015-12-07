#include "vertToCell.h"
#include <list>
#include <iterator>
#include <iostream>
#include "cellToVert.h"

using namespace std;

void vertToCell::push(int i){
	data.push_back(i);
}

void vertToCell::size(){
	data.size();
}

void vertToCell::show(){
	list<int>::iterator it = data.begin();
	while(it != data.end()){
		cout << *it << " ";
		it++;
	}
	cout << "\n";
}

int& vertToCell::operator[](int i){
	
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

	cout << "error 2075008437q-" << endl;
	
	return data.front();
}

void vertToCell::setDataBase(cellToVert* CV, int Ncell, int n){
	data.clear();
	for(int i=0; i<Ncell; i++){
		for(int j=0; j<(int)CV[i].data.size(); j++){
			if(n == CV[i][j]){
				data.push_back(i);
			}
		}
	}
}
