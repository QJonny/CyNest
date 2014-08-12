#include "tools.h"

// vector 3d
Vector3d::Vector3d(double x__, double y__, double z__) {
	x_ = x__;
	y_ = y__;
	z_ = z__;
	norm_ = sqrt(x_*x_ + y_*y_ + z_*z_);
}

Vector3d::Vector3d() {
	Vector3d(0.0, 0.0, 0.0);
}

void Vector3d::set(double x__, double y__, double z__) {
	x_ = x__;
	y_ = y__;
	z_ = z__;
	norm_ = sqrt(x_*x_ + y_*y_ + z_*z_);
}

double Vector3d::norm() {
	return norm_;
}

Vector3d Vector3d::add(Vector3d v) {
	return Vector3d(x_ + v.x_, y_ + v.y_, z_ + v.z_);
}

Vector3d Vector3d::sub(Vector3d v) {
	return Vector3d(x_ - v.x_, y_ - v.y_, z_ - v.z_);
}

Vector3d Vector3d::mul(double d) {
	return Vector3d(x_ * d, y_ * d, z_ * d);
}

Vector3d Vector3d::div(double d) {
	return Vector3d(x_ / d, y_ / d, z_ / d);
}

Vector3d Vector3d::normalize() {
	double n = norm();

	return Vector3d(x_ / n, y_ / n, z_ / n);
}

double Vector3d::x() {
	return x_;
}

double Vector3d::y() {
	return y_;
}

double Vector3d::z() {
	return z_;
}






// takes a string of the format "[val1,val2,val3,...]"
// and returns a double array of these values
bool parseList(char* str, double* list) {
	char valBuffer[LIST_ELEMENT_SIZE] = {0};
	char* cChar = &str[1]; // current char pointer, neglect the initial '['
	int j = 0; // valBuffer position
	int i = 0; // i: current value (val1, val2, ...)
	
	while(true) { 
		if(*cChar == ' ') {
			// nothing, skip
		}
		else if(*cChar == ',') {
			valBuffer[j] = '\0';
			j = 0;
			list[i] = atof(valBuffer);
			i++;
		}
		else if(*cChar == ']') {
			valBuffer[j] = '\0';
			list[i] = atof(valBuffer);
			break;
		}
		else if(*cChar >= 48 && *cChar <= 57 || *cChar == '-' || *cChar == '.') { // 0 - 9, -, .
			valBuffer[j] = *cChar;
			j++;
		}
		else {
			return false;
		}
		cChar++;
	}
	
	return true;
}



double generateRandomNumber(double low, double high) {
	return (high-low)*((float)rand()/RAND_MAX) + low;
}




string deleteExecName(string path) {
        int slashPointer = path.length() - 1;
        while(slashPointer >= 0 && path[slashPointer] != '/') {
                slashPointer--;
        }
        return path.substr(0, slashPointer + 1);
}

string getExecDirectory() {
  char pBuf[FILENAME_MAX];
#ifdef WINDOWS
        int bytes = GetModuleFileName(NULL, pBuf, FILENAME_MAX);
        if(bytes == 0)
                return NULL;
        else
                return deleteExecName(string(pBuf));

#else
        char szTmp[32];
        sprintf(szTmp, "/proc/%d/exe", getpid());
        int bytes = min((int)readlink(szTmp, pBuf, FILENAME_MAX), FILENAME_MAX - 1);
        if(bytes >= 0)
                pBuf[bytes] = '\0';
        return deleteExecName(string(pBuf));
#endif
}

