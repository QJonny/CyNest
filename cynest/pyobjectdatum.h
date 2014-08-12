/*
 *  pydatum.h
 *
 *  This file is part of NEST
 *
 *  Copyright (C) 2012 by
 *  The NEST Initiative
 *
 *  See the file AUTHORS for details.
 *
 *  Permission is granted to compile and modify
 *  this file for non-commercial use.
 *  See the file LICENSE for details.
 *
 */

#ifndef PYOBJECTDATUM_H
#define PYOBJECTDATUM_H

#include <Python.h>
#include "datum.h"

#include <string>
#include "dataconverter.h"

#define GET_STATUS_METHOD 1
#define SET_STATUS_METHOD 2

// prefixed all references to members of GenericDatum with this->,
// since HP's aCC otherwise complains about them not being declared
// according to ISO Standard Sec. 14.6.2(3) [temp.dep]
// HEP, 2001-08-08



class PyObjectDatum: public Datum
{

private:
  PyObject* pyObj;
  double* currents;
  double* in_spikes;
  double* ex_spikes;
  long* t_lag;
  long* spike;
  double* current_value;
  
  int forbiddenParamsLength;
  std::string forbiddenParams[16];

  PyCFunction updateFct;



struct PyMethodDef* getUpdateRef(struct PyMethodDef *tp_methods) {
	int i = 0;
	while(&(tp_methods[i]) != NULL) {
		if(std::string("update").compare(tp_methods[i].ml_name) == 0) {
			return &(tp_methods[i]);
		}
		i++;
	}
	return NULL;
}


bool isOK(std::string method)
{	
	for(int i = 0; i < forbiddenParamsLength; i++) {
		if(method.compare(forbiddenParams[i]) == 0) {
			return false;
		}
	}
	
	return true;
}

void init() {
	  currents = NULL;
	  in_spikes = NULL;
	  ex_spikes = NULL;
	  t_lag = NULL;
	  spike = NULL;
	  // the cython model shouldn't use or access these parameters
	  forbiddenParamsLength = 16;
	  forbiddenParams[0] = "archiver_length";
	  forbiddenParams[1] = "frozen";
	  forbiddenParams[2] = "global_id";
	  forbiddenParams[3] = "local";
	  forbiddenParams[4] = "local_id";
	  forbiddenParams[5] = "model";
	  forbiddenParams[6] = "node_type";
	  forbiddenParams[7] = "parent";
	  forbiddenParams[8] = "pyobject";
	  forbiddenParams[9] = "recordables";
	  forbiddenParams[10] = "state";
	  forbiddenParams[11] = "t_spike";
	  forbiddenParams[12] = "tau_minus";
	  forbiddenParams[13] = "tau_minus_triplet";
	  forbiddenParams[14] = "thread";
	  forbiddenParams[15] = "vp";
}

public:
  PyObjectDatum(PyObject *py) {
	  pyObj = py;

	  Py_XINCREF(pyObj);

	  init();
  }

  PyObjectDatum(PyObjectDatum* py) {
	  this->pyObj = py->pyObj;

	init();
  }
  ~PyObjectDatum() {}

  PyObjectDatum *clone(void) const
  {
        return new PyObjectDatum(*this);
  }

void putStdParams(double** curr, double** is, double** es, long** tl, long** sp, double** cv) {
	// important, otherwise segmentation fault
	PyGILState_STATE s = PyGILState_Ensure();
	
	// numeric conversion in order to create a long variable containing the address of the std params 
	*curr = currents = PyInt_AsLong(PyObject_CallMethod(this->pyObj, "getPCurrents", NULL));
	*is = in_spikes = PyInt_AsLong(PyObject_CallMethod(this->pyObj, "getPIn_Spikes", NULL));
	*es = ex_spikes = PyInt_AsLong(PyObject_CallMethod(this->pyObj, "getPEx_Spikes", NULL));
	*tl = t_lag = PyInt_AsLong(PyObject_CallMethod(this->pyObj, "getPT_Lag", NULL));
	*sp = spike = PyInt_AsLong(PyObject_CallMethod(this->pyObj, "getPSpike", NULL));
	*cv = current_value = PyInt_AsLong(PyObject_CallMethod(this->pyObj, "getPCurrent_Value", NULL));
	
	// we also extract the direct pointer to the update c function
	struct PyMethodDef* updateRef = getUpdateRef(this->pyObj->ob_type->tp_methods);
	if(updateRef != NULL) {
		this->updateFct = updateRef->ml_meth;
	}

	PyGILState_Release(s);
}

void call_method(std::string cmd) {
	// important, otherwise segmentation fault
	PyGILState_STATE s = PyGILState_Ensure();

	PyObject_CallMethod(this->pyObj, cmd.c_str(), NULL);

    PyGILState_Release(s);
}

void call_update() {
	// important, otherwise segmentation fault
	PyGILState_STATE s = PyGILState_Ensure();
	
	updateFct(this->pyObj, NULL);

    PyGILState_Release(s);
}

void call_update_optimized() {
	// without the GIL the function is faster and, much more
	// important, the multithreading is not affected
	updateFct(this->pyObj, NULL);
}



void call_status_method(int m, void* status_) {
	// important, otherwise segmentation fault
	PyGILState_STATE s = PyGILState_Ensure();
	DictionaryDatum* status = static_cast<DictionaryDatum*>(status_);
	
	if(m == SET_STATUS_METHOD) {
		// creation of an empty python dictionary
		PyObject* dict = PyDict_New();

		// filling the dictionary
		for(Dictionary::iterator it = (*status)->begin(); it != (*status)->end(); ++it) {
			// if the status_ element is not forbidden (is actually one of the model parameters), it is copied
			if(isOK(it->first.toString()) == true) {
			#if PY_MAJOR_VERSION >= 3
				PyObject_SetItem(dict, PyUnicode_FromString(it->first.toString().c_str()), dataConverter.datumToObject( (**status)[it->first.toString()].datum()));
			#else
				PyObject_SetItem(dict, PyString_FromString(it->first.toString().c_str()), dataConverter.datumToObject( (**status)[it->first.toString()].datum()));				
			#endif
			}
		}
		
		// after the python dict has been filled, we can call the method
		#if PY_MAJOR_VERSION >= 3
			PyObject_CallMethodObjArgs(this->pyObj, PyUnicode_FromString("setStatus"), dict, NULL);
		#else
			PyObject_CallMethodObjArgs(this->pyObj, PyString_FromString("setStatus"), dict, NULL);
		#endif
                Py_XDECREF(dict);
	}
	else if(m == GET_STATUS_METHOD) {
		PyObject* dict = PyObject_CallMethod(this->pyObj, "getStatus", NULL);
		PyObject* subPyObj=0;
		PyObject *key=0;
		Py_ssize_t pos = 0;

		// looping through the received python dictionary
		while (PyDict_Next(dict, &pos, &key, &subPyObj)) 
		{
		#if PY_MAJOR_VERSION >= 3
			PyObject* pStrObj = PyUnicode_AsUTF8String(key); 
			char* zStr = PyBytes_AsString(pStrObj); 
			Py_DECREF(pStrObj);
			(**status)[zStr] = dataConverter.objectToDatum(subPyObj);
		#else
			(**status)[PyString_AsString(key)] = dataConverter.objectToDatum(subPyObj);
		#endif
		}

		Py_XDECREF(dict);
		// the std params also have to be updated
		if(currents != NULL && in_spikes != NULL && ex_spikes != NULL && t_lag != NULL && spike != NULL) {
			(**status)["currents"] = *currents;
			(**status)["in_spikes"] = *in_spikes;
			(**status)["ex_spikes"] = *ex_spikes;
			(**status)["t_lag"] = *t_lag;
			(**status)["spike"] = *spike;
			(**status)["current_value"] = *current_value;
		}
	}
	PyGILState_Release(s);
}

  /**
   * Accept a DatumVisitor as a visitor to the datum (visitor pattern).
   * This member has to be overridden in the derived classes
   * to call visit and passing themselves as an argument.
   */
  void use_converter(DatumConverter &v) {  }
  
  void  print(std::ostream & ) const {printf("This is a PyObjectDatum");}
  void  pprint(std::ostream &) const {printf("This is a PyObjectDatum");}

  void  list(std::ostream &o, std::string prefix, int l) const
  {
    if(l==0)
      prefix="-->"+prefix;
    else
      prefix="   "+prefix;
    o << prefix;
    print(o);
  }

  void  input_form(std::ostream &o) const
  {
    pprint(o);
  }
  
  bool  equals(const Datum *d) const
  {
    return this == d;
  }
  
  void  info(std::ostream &) const{printf("PyObjectDatum");}
};



#endif
