#include "vertex.h"
#include <iomanip>

vertex::vertex(){
	x = 0.0;
	y = 0.0;
	z = 0.0;
	b = 0;
}

vertex::vertex(double a, double b){
	x = a;
	y = b;
	z = 0.0;
}

vertex::vertex(double a, double b, double c){
	x = a;
	y = b;
	z = c;
}

void vertex::set(double a, double b){
	x = a;
	y = b;
}

vertex vertex::operator+(vertex obj){
	vertex ans;
	ans.x = this->x + obj.x;
	ans.y = this->y + obj.y;
	ans.z = this->z + obj.z;
	return ans;
}

vertex vertex::operator-(vertex obj){
	vertex ans;
	ans.x = this->x - obj.x;
	ans.y = this->y - obj.y;
	ans.z = this->z - obj.z;
	return ans;
}

vertex vertex::operator-(){
	vertex ans;
	ans.x = - this->x;
	ans.y = - this->y;
	ans.z = - this->z;
	return ans;
}

vertex vertex::operator/(double a){
	vertex ans;
	ans.x = this->x/a;
	ans.y = this->y/a;
	ans.z = this->z/a;
	return ans;
}
	
double vertex::operator*(vertex obj){
	return this->x * obj.x + this->y * obj.y + this->z * obj.z;
}

vertex vertex::operator%(vertex obj){
	vertex ans;
	ans.x = y * obj.z - z * obj.y;
	ans.y = z * obj.x - x * obj.z;
	ans.z = x * obj.y - y * obj.x;
	return ans;
}

vertex operator*(double a, vertex obj){
	vertex ans;
	ans.x = a * obj.x;
	ans.y = a * obj.y;
	ans.z = a * obj.z;
	return ans;
}

void vertex::show(){
	cout << "x:" << setprecision(14) << x << endl;
	cout << "y:" << setprecision(14) << y << endl;
	cout << "z:" << setprecision(14) << z << endl;
}

double vertex::norm(){
	return sqrt(x*x + y*y + z*z);
}

double th = M_PI / 2.0;

void vertex::rot(double r){

	double temp_x = 0;
	double temp_y = 0;

	temp_x = cos(th) * x - sin(th) * y;
	temp_y = sin(th) * x + cos(th) * y;

	x = r * temp_x;
	y = r * temp_y;
}

/*-------------------------------------------*/

vertex cross(vertex v1, vertex v2){
	vertex ans;
	ans.x = v1.y*v2.z - v1.z*v2.y;
	ans.y = v1.z*v2.x - v1.x*v2.z;
	ans.z = v1.x*v2.y - v1.y*v2.x;
	return ans;
}

vertex av(double a, vertex v){
	vertex ans;
	ans.x = a * v.x;
	ans.y = a * v.y;
	ans.z = a * v.z;
	return ans;
}
