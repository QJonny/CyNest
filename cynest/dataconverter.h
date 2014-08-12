#ifndef DATACONVERTER_H
#define DATACONVERTER_H

#include <Python.h>
#include "cynestkernel.h"
#include "datumtopythonconverter.h"

class DataConverter
{
private:
    DatumToPythonConverter dTp;
    NESTEngine pTd;

public:
	DataConverter() {
		dTp = DatumToPythonConverter();
		pTd = NESTEngine();
	}
    PyObject* datumToObject(Datum* d){
        return dTp.convertDatum(d);
	}

    Datum* objectToDatum(PyObject* o){
        return pTd.PyObject_as_Datum(o);
	}
};

DataConverter dataConverter = DataConverter();

#endif
