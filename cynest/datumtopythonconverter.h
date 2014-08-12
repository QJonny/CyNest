/*
 *  datumtopythonconverter.h
 *
 *  This file is part of NEST
 *
 *  Copyright (C) 2008 by
 *  The NEST Initiative
 *
 *  See the file AUTHORS for details.
 *
 *  Permission is granted to compile and modify
 *  this file for non-commercial use.
 *  See the file LICENSE for details.
 *
 */

#ifndef DATUMTOPYTHONCONVERTER_H
#define DATUMTOPYTHONCONVERTER_H

#include<Python.h>
#include "py3k_compatibility.h"

#ifdef HAVE_NUMPY

// These defines are needed for Python C extensions using NumPy
//
// Since we have an extension consisting of several source files,
// we need to define PY_ARRAY_UNIQUE_SYMBOL to the same symbol
// in all files (according to the numpy documentation).
//
// Unless stated in the documentation, we also have to define
// NO_IMPORT_ARRAY in all files, which do not call import_array.
// import_array is only called in pynestkernel.cpp.
// (this becomes clear from reading __multiarray_api.h). Omitting
// this will cause a segfault.
// 
// Numpy uses the unique symbol do define and initialize a pointer
// void *PY_ARRAY_UNIQUE_SYMBOL in the file calling import_array.
// All other files have an extern void* PY_ARRAY_UNIQUQ_SYMBOL
// definition.
//
// Moritz, 2007-04-17

#define NO_IMPORT_ARRAY
#define PY_ARRAY_UNIQUE_SYMBOL _pynest_arrayu

extern "C" {
#include <numpy/arrayobject.h>
}

#endif

#include "dictdatum.h"
#include "dictutils.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "arraydatum.h"
#include "booldatum.h"
#include "stringdatum.h"
#include "namedatum.h"
#include "connectiondatum.h"
#include "datumconverter.h"

/**
 * Converter class for SLI Datums to python datums.
 * This is a DatumVisitor (visitor pattern), which visits
 * the SLI Datum. Each SLI Datum will call visit(*this) and therefore evokes
 * the approproate (overloaded) implementation of visit here.
 */
class DatumToPythonConverter : public DatumConverter
{
 public:
  DatumToPythonConverter();
  
  /**
   * For every datum type to be converted we have to override the default visit implementation
   * of DatumVisitor with the respective type of the datum to be converted.
   */
  void convert_me(Datum &d);
  void convert_me(DoubleDatum &d);
  void convert_me(IntegerDatum &i);
  void convert_me(BoolDatum &i);
  void convert_me(StringDatum &s);
  void convert_me(ArrayDatum &ad); 
  void convert_me(DictionaryDatum &dd);
  void convert_me(LiteralDatum &ld);
  void convert_me(DoubleVectorDatum &dvd);
  void convert_me(IntVectorDatum &dvd);
  void convert_me(ConnectionDatum &dvd);

  /**
   * Converts a given NEST Datum to a Python object. The return value is a new reference to a Python object.
   * The user is responsible for disposing the reference at the end of the program.
   */
  PyObject * convert(Datum &d);

  // Added by Jonny Quarta
  PyObject * convertDatum(Datum *d);

  inline
  PyObject *getPyObject()
  {
    return py_object_;
  }

 private:
  PyObject *py_object_;

};

#endif
