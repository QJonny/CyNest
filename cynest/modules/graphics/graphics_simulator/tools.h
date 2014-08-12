#ifndef TOOLS_H
#define TOOLS_H


#include "headers.h"


using namespace std;


class Vector3d {
private:
	double x_;
	double y_;
	double z_;
	double norm_;

public:	
	Vector3d();
	
	Vector3d(double x__, double y__, double z__);
	
	void set(double x__, double y__, double z__);
	
	double norm();
	Vector3d add(Vector3d v);
	Vector3d sub(Vector3d v);
	Vector3d mul(double d);
	Vector3d div(double d);
	Vector3d normalize();
	double x();
	double y();
	double z();
};




bool parseList(char* str, double* list);

double generateRandomNumber(double low, double high);

string getExecDirectory();
string deleteExecName(string path);

#endif


