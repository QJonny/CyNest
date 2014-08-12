#ifndef OBJECT_MANAGER_H
#define OBJECT_MANAGER_H

#include "../nestkernel/nest_time.h"
#include "../nestkernel/scheduler.h"

using namespace nest;

class UnitManager {
private:
	/*
	 1 tic
	 2 step
	 3 ms
	 4 ms_stamp
	*/
	int unit;
	long lValue;
	double dValue;

public:
	UnitManager(int u, double value) {
		this->unit = u;

		this->dValue = value;
	}

	UnitManager(int u, long value) {
		this->unit = u;

		this->lValue = value;
	}

	Time generateTime() {
		if(this->unit == 1) {
			return Time(Time::tic(this->lValue));
		} else if(this->unit == 2) {
			return Time(Time::step(this->lValue));
		} else if(this->unit == 3) {
			return Time(Time::ms(this->dValue));
		} else if(this->unit == 4) {
			return Time(Time::ms_stamp(this->dValue));
		}
		return Time(Time::ms(0.0));
	}
};

#endif
