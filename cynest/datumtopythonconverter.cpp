
/*
 *  datumtopythonconverter.h
 *
 *  This file is part of NEST
 *
 *  Copyright (C) 2004 by
 *  The NEST Initiative
 *
 *  See the file AUTHORS for details.
 *
 *  Permission is granted to compile and modify
 *  this file for non-commercial use.
 *  See the file LICENSE for details.
 *
 *
 */

#include "datumtopythonconverter.h"

#include "kernel.h"
void DatumToPythonConverter::convert_me(Datum &s)
{
    Token *t= new Token(s); 
    py_object_ = Token_to_PyObject(t);
}

DatumToPythonConverter::DatumToPythonConverter()
{
  py_object_ = 0;
}

void DatumToPythonConverter::convert_me(DoubleDatum &d)
{
  py_object_ = PyFloat_FromDouble(d.get());
}

void DatumToPythonConverter::convert_me(IntegerDatum &i)
{
  py_object_ = PyInt_FromLong(i.get());
}
  
void DatumToPythonConverter::convert_me(BoolDatum &b)
{
  
  if (b.get()) { // true
    Py_INCREF(Py_True);
    py_object_ = Py_True;
  }
  else { // false
    Py_INCREF(Py_False);
    py_object_ = Py_False;
  }

}

void DatumToPythonConverter::convert_me(StringDatum &s)
{
#if PY_MAJOR_VERSION >= 3
    const char error[]="Can't encode SLI string to PyUnicode";
    py_object_ = PyUnicode_DecodeUTF8(s.c_str(),s.size(), error);
#else
    py_object_ = PyString_FromString(s.c_str());
#endif
}


void DatumToPythonConverter::convert_me(ArrayDatum &ad)
{
  py_object_ = PyList_New(ad.size());

  DatumToPythonConverter dpc;

  int i = 0;
  for(Token *idx = ad.begin(); idx != ad.end(); ++idx) {
    // recurse to convert this element in the array
    // add it to the python list
    PyList_SetItem(py_object_, i, dpc.convert(*idx->datum()));
    ++i;
  }
}

void DatumToPythonConverter::convert_me(DictionaryDatum &dd)
{  
    py_object_ = PyDict_New();    
    const Token* subt;

    DatumToPythonConverter dpc;

    for(TokenMap::const_iterator where = dd->begin(); where != dd->end(); where ++) {
      subt = (Token*) ( &((*where).second) );

      PyDict_SetItemString(py_object_, (*where).first.toString().c_str(), dpc.convert(*subt->datum()));
      // The reference in dcp is no longer needed, so we dispose it.
      Py_DECREF(dpc.getPyObject());
    }
}

void DatumToPythonConverter::convert_me(LiteralDatum &ld)
{
    std::string s=ld.toString();
#if PY_MAJOR_VERSION >= 3
    const char error[]="Can't encode SLI literal to PyUnicode";
    py_object_ = PyUnicode_DecodeUTF8(s.c_str(),s.size(), error);
#else
    py_object_ = PyString_FromString(s.c_str());
#endif
}

void DatumToPythonConverter::convert_me(DoubleVectorDatum &dvd)
{
  int dims = dvd->size();

#ifdef HAVE_NUMPY
  PyArrayObject *array;

// PyArray_SimpleNew is a drop-in replacement for PyArray_FromDims
#if (NPY_VERSION >= 0x01000009)
// PyArray_SimpleNew takes ``npy_intp*`` dims instead of ``int*`` dims
// which matters on 64-bit systems and it does not initialize the
// memory to zero.
  npy_intp npydims = dims;
  array = (PyArrayObject*)(PyArray_SimpleNew(1, &npydims, PyArray_DOUBLE));
#else
  array = (PyArrayObject*)(PyArray_FromDims(1, &dims, PyArray_DOUBLE));  
#endif 

  std::copy( dvd->begin(), dvd->end(), reinterpret_cast<double*>(array->data) );  
  py_object_ = (PyObject*) array;
#else
  py_object_ = PyList_New(dims);
  for(int i=0; i<dims; i++)
    PyList_SetItem(py_object_, i, PyFloat_FromDouble(dvd->at(i)));
#endif //HAVE_NUMPY
}

void DatumToPythonConverter::convert_me(IntVectorDatum &ivd)
{
  int dims = ivd->size();

#ifdef HAVE_NUMPY
  PyArrayObject *array;

// PyArray_SimpleNew is a drop-in replacement for PyArray_FromDims
#if (NPY_VERSION >= 0x01000009)
// PyArray_SimpleNew takes ``npy_intp*`` dims instead of ``int*`` dims
// which matters on 64-bit systems and it does not initialize the
// memory to zero.
  npy_intp npydims = dims;
  array = (PyArrayObject*)PyArray_SimpleNew(1, &npydims, PyArray_INT);
#else
  array = (PyArrayObject*)PyArray_FromDims(1, &dims, PyArray_INT);
#endif

  std::copy( ivd->begin(), ivd->end(), reinterpret_cast<int*>(array->data) );
  py_object_ = (PyObject*) array;
#else
  py_object_ = PyList_New(dims);
  for(int i=0; i<dims; i++)
    PyList_SetItem(py_object_, i, PyInt_FromLong(ivd->at(i)));
#endif //HAVE_NUMPY
    
}

void DatumToPythonConverter::convert_me(ConnectionDatum &cd)
{
    const int dims=5;
#ifdef HAVE_NUMPY
  PyArrayObject *array;

// PyArray_SimpleNew is a drop-in replacement for PyArray_FromDims
#if (NPY_VERSION >= 0x01000009)
// PyArray_SimpleNew takes ``npy_intp*`` dims instead of ``long*`` dims
// which matters on 64-bit systems and it does not initialize the
// memory to zero.
npy_intp npydims = dims;

  array = (PyArrayObject*)PyArray_SimpleNew(1, &npydims, PyArray_LONG);
#else
  array = (PyArrayObject*)PyArray_FromDims(1, &dims, PyArray_LONG);
#endif
  long *val=reinterpret_cast<long*>(array->data);
  *val++= cd.get_source_gid();
  *val++= cd.get_target_gid();
  *val++= cd.get_target_thread();
  *val++= cd.get_synapse_model_id();
  *val++= cd.get_port();

  py_object_ = (PyObject*) array;
#else
  py_object_ = PyList_New(dims);
  int i=0;
  PyList_SetItem(py_object_, i++, PyInt_FromLong(cd.get_source_gid()));
  PyList_SetItem(py_object_, i++, PyInt_FromLong(cd.get_target_gid()));
  PyList_SetItem(py_object_, i++, PyInt_FromLong(cd.get_target_thread()));
  PyList_SetItem(py_object_, i++, PyInt_FromLong(cd.get_synapse_model_id()));
  PyList_SetItem(py_object_, i++, PyInt_FromLong(cd.get_port()));
#endif //HAVE_NUMPY
}

PyObject *DatumToPythonConverter::convert(Datum &d)
{
  d.use_converter(*this);
  return py_object_;
}

// Added by Jonny Quarta
PyObject *DatumToPythonConverter::convertDatum(Datum *d)
{
  return convert(*d);
}



