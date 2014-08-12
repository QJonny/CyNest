/*
 *  cynestkernel.cpp
 *
 *  Interface between Python and the NEST simulation tool. 
 *  www.nest-initiative.org
 *
 *  Developed in the context of the EU FACETS project.
 *  facets.kip.uni-heidelberg.de
 *
 *  Authors:
 *  Marc-Oliver Gewaltig
 *  Eilif Muller
 *  Moritz Helias
 *  Jochen Martin Eppler
 *
 *  This file is part of NEST
 *
 *  Copyright (C) 2004-2009 by
 *  The NEST Initiative
 *
 *  See the file AUTHORS for details.
 *
 *  Permission is granted to compile and modify
 *  this file for non-commercial use.
 *  See the file LICENSE for details.
 */


#include <Python.h>

// If we're not using python 2.5
#if (PY_VERSION_HEX < 0x02050000)
typedef int Py_ssize_t;
#endif

// Needed for Python C extensions using NumPy, see comment in
// datumtopythonconverter.h
#define PY_ARRAY_UNIQUE_SYMBOL _pynest_arrayu

#ifdef HAVE_NUMPY
extern "C"
{
#include <numpy/arrayobject.h>
}
#endif

#include "interpret.h"
#include "network.h"
#include "communicator.h"
#include "random_numbers.h"
#include "slistartup.h"
#include "sliarray.h"
#include "processes.h"
#include "nestmodule.h"
#include "dynamicloader.h"
#include "oosupport.h"
#include "processes.h"
#include "sliregexp.h"
#include "specialfunctionsmodule.h"
#include "sligraphics.h"
#include "compose.hpp"
#include "filesystem.h"

#include "doubledatum.h"
#include "integerdatum.h"
#include "dictdatum.h"
#include "dictutils.h"
#include "arraydatum.h"
#include "booldatum.h"
#include "stringdatum.h"
#include "pyobjectdatum.h"
#include "connectiondatum.h"
#include "datumtopythonconverter.h"
#include "psignal.h"

#include <algorithm>

#include "spikecounter.h"

#include "static_modules.h"

#include "cynestkernel.h"
extern "C" {
#include "kernel.h"
}

// Little hack for making the file recognized. Otherwise,
// the compiler doesn't seem to consider it (only the .h) 
#include "cython_neuron.cpp"

/*
The following instance of nest::spikecounter needs to be defined
under MacOS X to prevent the linker from throwing away the
constructor of spikecounter during linking.
See also bug #301.
MH, 2009/01/07
Update Nov. 2012
At least on Mac OS 10.7 (lion) this declaration is no longer needed.
I commented it our. If problems don't re-occur, we can remove the code
in future versions. mog
*/
#ifdef __APPLE__
//http://developer.apple.com/technotes/tn2002/tn2071.html#Section10
//nest::spikecounter pseudo_spikecounter_instance(0.0,0.0);
#endif

extern "C"
{
void init_numpy_()
#ifdef HAVE_NUMPY
//  we need to set up the numeric array type
  import_array()
#else
{}
#endif
}

NESTEngine::NESTEngine()
    : initialized_(false),
      NESTError_(0), //!< Python error object. We are responsible.
      pEngine_(0),   //!< NEST instance. We are responsible for this pointer
      pNet_(0)       //!< Pointer given to NEST. NEST will delete this object
{

}

NESTEngine::~NESTEngine()
{
    if(initialized_)
    {
	delete pEngine_;
	Py_XDECREF(NESTError_);
    }
}

/**
 * Initialize numpy if it is available. 
 * This is a stand-alone function, because import_array() seems to be a macro that resolves to a return 
 * statement. 
 */
void NESTEngine::init_numpy()
{ // we need to set up the numeric array type
    init_numpy_();
}


bool NESTEngine::check_engine()
{
    if (not initialized_ or pEngine_==NULL) 
    {
	PyErr_SetString(NESTError_, "runsli(): PyNEST engine not initialized properly or finalized already.");
	return false;
    }
    return true;
}

bool NESTEngine::init(std::vector<string> argv, std::string modulepath)
{
    int argc=argv.size();
    // We create a Python exception object. This object is used by all functions.
    // Calling PyErr_SetString will trigger the exception.
    static char nest_error_name[]="cynest.NESTError";
    NESTError_ = PyErr_NewException(nest_error_name,NULL,NULL);

    init_numpy();
    
    if (initialized_)
    {
	PyErr_SetString(NESTError_, "NESTEngine::init(): NEST is already initialized.");
	return false;
    }
    
    pEngine_ = new SLIInterpreter;
    
    if (pEngine_==NULL)
    {
	PyErr_SetString(NESTError_, "NESTEngine::init(): Cannot create NEST interpreter!");
	return false;
    }
    
    char **c_argv=new char *[argc];
    for(int i=0; i< argc; ++i)
	c_argv[i] = const_cast<char *>(argv[i].c_str());

#ifdef HAVE_MPI
    nest::Communicator::init(&argc, &c_argv);
#endif
    
    addmodule<SLIArrayModule>(*pEngine_);
    addmodule<OOSupportModule>(*pEngine_);
    addmodule<RandomNumbers>(*pEngine_);
    addmodule<SpecialFunctionsModule>(*pEngine_);   // safe without GSL
    addmodule<SLIgraphics>(*pEngine_);
    pEngine_->addmodule(new SLIStartup(argc,c_argv));
    addmodule<Processes>(*pEngine_);
    addmodule<RegexpModule>(*pEngine_);
    addmodule<FilesystemModule>(*pEngine_);

  // create the network and register with NestModule class
    pNet_ = new nest::Network(*pEngine_);
    assert(pNet_ != 0);
    nest::NestModule::register_network(*pNet_);
    addmodule<nest::NestModule>(*pEngine_);

  // now add static modules providing models
    add_static_modules(*pEngine_, *pNet_);

#ifdef HAVE_LIBLTDL  // no dynamic loading without LIBLTDL
    //dynamic loader module for managing linked and dynamically loaded extension modules
    nest::DynamicLoaderModule *pDynLoader = new nest::DynamicLoaderModule(pNet_, *pEngine_);

  // initialize all modules that were linked into at compile time
  // these modules have registered via calling DynamicLoader::registerLinkedModule
  // from their constructor
  pDynLoader->initLinkedModules(*pEngine_);

  // interpreter will delete module on destruction
  pEngine_->addmodule(pDynLoader);
#endif

  // add the init-script to the list of module initializers
  ArrayDatum *ad = dynamic_cast<ArrayDatum *>(pEngine_->baselookup(pEngine_->commandstring_name).datum());
  assert(ad != NULL);
  ad->push_back(new StringDatum("(" + modulepath + "/pynest-init.sli) run"));

  // clean up argv memory
  delete [] c_argv;
  pEngine_->startup();
  initialized_=true;

 
  return initialized_;
}

bool NESTEngine::push(PyObject *args)
{
    if( not check_engine())
	return false;

    Datum *pdat = PyObject_as_Datum(args);
    if (pdat != 0) {
	pEngine_->OStack.push(pdat);
	return true;
    }
    else
	return false;
}

bool NESTEngine::push_token(Token obj)
{
    if( not check_engine() or  obj.empty())
	return false;

    pEngine_->OStack.push_move(obj);
    return true;
}

bool NESTEngine::push_connections(PyObject *arg)
{
    if(not check_engine())
	return false;

    if (!PyList_Check(arg) && !PyTuple_Check(arg))
    {
	PyErr_SetString(NESTError_, "NESTEngine::push_connections(): Argument must be a list of dictionaries.");
	return false;
    }

    ArrayDatum connectome;
    PyObject* subPyObj;

    size_t size = (PyList_Check(arg)) ? PyList_Size(arg) : PyTuple_Size(arg);
    connectome.reserve(size);
 
    for (size_t i = 0; i < size; ++i)
    {
	subPyObj = (PyList_Check(arg)) ? PyList_GetItem(arg, i) : PyTuple_GetItem(arg, i);
	
	if (PyDict_Check(subPyObj))
	{
	    PyObject* subsubPyObj=PyDict_GetItemString(subPyObj, nest::names::source.toString().c_str());
	    long source;
	    long target_thread;
	    long synapse_modelid;
	    long port;
	    
	    if (subsubPyObj != NULL && PyInt_Check(subsubPyObj))
		source = PyInt_AsLong(subsubPyObj);
	    else
	    {
		PyErr_SetString(NESTError_, "NESTEngine::push_connections(): No source entry in dictionary.");
		return false;
	    }
	    
	    subsubPyObj = PyDict_GetItemString(subPyObj, nest::names::target_thread.toString().c_str());
	    if (subsubPyObj != NULL && PyInt_Check(subsubPyObj))
		target_thread = PyInt_AsLong(subsubPyObj);
	    else
	    {
		PyErr_SetString(NESTError_, "NESTEngine::push_connections(): No target_thread entry in dictionary.");
		return false;
	    }
	    
	    subsubPyObj = PyDict_GetItemString(subPyObj, nest::names::synapse_modelid.toString().c_str());
	    if (subsubPyObj != NULL && PyInt_Check(subsubPyObj))
		synapse_modelid = PyInt_AsLong(subsubPyObj);
	    else
	    {
		PyErr_SetString(NESTError_, "NESTEngine::push_connections(): No synapse_modelid entry in dictionary.");
		return false;
	    }
	    
	    subsubPyObj = PyDict_GetItemString(subPyObj, nest::names::port.toString().c_str());
	    if (subsubPyObj != NULL && PyInt_Check(subsubPyObj))
		port = PyInt_AsLong(subsubPyObj);
	    else
	    {
		PyErr_SetString(NESTError_, "NESTEngine::push_connections(): No port entry in dictionary.");
		return false;
	    }
	    
	    ConnectionDatum cd = ConnectionDatum(nest::ConnectionID(source, target_thread, synapse_modelid, port));
	    connectome.push_back(cd);
	}
#ifdef HAVE_NUMPY
	else if (PyArray_Check(subPyObj) )
	{
	    size_t array_size = PyArray_Size(subPyObj);
	    if(array_size!=5)
	    {
		std::string error = String::compose("push_connections(): At position %1 in connection ID list.",i)+
		    "\n Connection ID must have exactly five entries.";
		PyErr_SetString(NESTError_, error.c_str());
		return false;
	    }
	    PyArrayObject *array = (PyArrayObject*) subPyObj;
	    assert(array != 0);
	    switch (array->descr->type_num)
	    {
	    case NPY_INT :
	    {
		PyArrayIterObject *iter= (PyArrayIterObject *)PyArray_IterNew(subPyObj);
		assert(iter != 0);

		int con[5];
		while(iter->index < iter->size)
		{
		    con[iter->index]= *(int*)(iter->dataptr);
		    PyArray_ITER_NEXT(iter);
		}
		delete iter ;
		connectome.push_back(new ConnectionDatum(nest::ConnectionID(con[0],con[1],con[2], con[3], con[4])));
		continue;
	    }
	    case NPY_LONG:
	    {
		PyArrayIterObject *iter= (PyArrayIterObject *)PyArray_IterNew(subPyObj);
		assert(iter != 0);

		long con[5];
		while(iter->index < iter->size)
		{
		    con[iter->index]= *(long*)(iter->dataptr);
		    PyArray_ITER_NEXT(iter);
		}
		delete iter ;
		connectome.push_back(new ConnectionDatum(nest::ConnectionID(con[0],con[1],con[2], con[3], con[4])));
		continue;
	    }
	    default:
		std::string error = String::compose("push_connection_arrays(): At position %1 in connection ID list.",i)+
		    "\n Connection ID must be a list or numpy array of five integers.";
		PyErr_SetString(NESTError_, error.c_str());
		return false;
	    }
	}
#endif
	else if (PyList_Check(subPyObj))
	{
	    size_t tuple_size = PyList_Size(subPyObj);
	    //std::cerr << tuple_size << std::endl;
	    if (tuple_size !=5)
	    {
		std::string error = String::compose("push_connection_arrays(): At position %1 in connection ID list.",i)+
		    "\n Connection ID must have exactly five entries.";
		PyErr_SetString(NESTError_, error.c_str());
		continue;
	    }
	    long con[5];
	    for (long j=0; j<5; ++j)
	    {
		PyObject* itemPyObj = PyList_GetItem(subPyObj, j);
		if(PyInt_Check(itemPyObj))
		{
		    con[j]=PyInt_AsLong(itemPyObj);
		}
#ifdef HAVE_NUMPY
		else if (PyArray_CheckScalar(itemPyObj)) // handle numpy array scalars
		{
		    PyArray_Descr *typecode;
		    typecode = PyArray_DescrFromScalar(itemPyObj);
		    switch (typecode->type_num)
		    {
		    case NPY_INT:
		    case NPY_LONG:
			long val;
			PyArray_ScalarAsCtype(itemPyObj, &val);
			con[j]=val;
			break;
		    default:
			std::string error = String::compose("push_connection_arrays(): At position %1: Unsupported Numpy array scalar type: '%2'.\n",
							    i,typecode->type_num);
			PyErr_SetString(NESTError_, error.c_str());
			return false;
		    }
		}
#endif
		else
		{
		    std::string error = String::compose("push_connection_arrays(): At position %1, %2 in connection ID list."
							"\n Connection ID must be a list, tuple, or and array of five integers.",i,j);
		    PyErr_SetString(NESTError_, error.c_str());
		    return false;
		}
	    }
	    connectome.push_back(new ConnectionDatum(nest::ConnectionID(con[0],con[1],con[2], con[3], con[4])));
	    continue;
	}
	else {
	    std::string error = String::compose("push_connection_arrays(): At position %1 in connection ID list.",i)+
		"\n Connection ID must be a list, tuple, or array of five integers.";
	    PyErr_SetString(NESTError_, error.c_str());
	    return false;
	}
    }
    
    pEngine_->OStack.push(connectome);
    return true;
}

bool NESTEngine::run(std::string cmd)
{
    if(not check_engine())
	return false;

    // send string to sli interpreter
    Py_BEGIN_ALLOW_THREADS
    // Safe the Python Signal handler
    Sigfunc *py_signal_handler= posix_signal(SIGINT, (Sigfunc *)SIG_IGN);
    // Set the SLISignal handler
    posix_signal(SIGINT,(Sigfunc *)SLISignalHandler);
    pEngine_->execute(cmd);
    // and now we re-set the Python signal handler
    posix_signal(SIGINT,(Sigfunc *)py_signal_handler);
    Py_END_ALLOW_THREADS
    
    return true;
}

bool NESTEngine::run_token(Token cmd)
{
    if(not check_engine() or not cmd)
	return false;

    // send string to sli interpreter
    Py_BEGIN_ALLOW_THREADS
    Sigfunc *py_signal_handler= posix_signal(SIGINT, (Sigfunc *)SIG_IGN);
    posix_signal(SIGINT,(Sigfunc *)SLISignalHandler);
    pEngine_->execute(cmd);
    posix_signal(SIGINT,(Sigfunc *)py_signal_handler);
    Py_END_ALLOW_THREADS
    
    return true;
}

PyObject *NESTEngine::pop()
{
    if(not check_engine())
	return false;

    if (pEngine_->OStack.empty()) 
    {
	PyErr_SetString(NESTError_, "NESTEngine::pop(): SLI stack is empty.");
	return NULL;
    }

    Token &t = pEngine_->OStack.top();
    DatumToPythonConverter DatumToPyObj;
    PyObject *pObj = 0;

    try {
	pObj = DatumToPyObj.convert(*t); // operator* returns reference to Datum
    }
    catch(TypeMismatch e)
    {
	PyErr_SetString(NESTError_, "NEST object cannot be converted to python object.");
    }
    pEngine_->OStack.pop();

    return pObj;
}

Token* NESTEngine::pop_token()
{
    if(not check_engine())
	return false;

    if (pEngine_->OStack.empty()) 
    {
	PyErr_SetString(NESTError_, "NESTEngine::pop(): SLI stack is empty.");
	return NULL;
    }

    Token *t= new Token();
    t->move(pEngine_->OStack.top());
    pEngine_->OStack.pop();
    
    return t;
}

void NESTEngine::register_cython_model(std::string model)
{
	nest::register_cython_model(pNet_, model);
}

Datum* NESTEngine::PyObject_as_Datum(PyObject *pObj)
{
  if (PyInt_Check(pObj)) { // object is integer or bool
    if (pObj==Py_True)
      return new BoolDatum(true);
    else if (pObj==Py_False)
      return new BoolDatum(false);
    else
      return new IntegerDatum(PyInt_AsLong(pObj));
  }

  if (PyFloat_Check(pObj)) // object is float
    return new DoubleDatum(PyFloat_AsDouble(pObj));

  if (PyString_Check(pObj)) // object is string
    return new StringDatum(PyString_AsString(pObj));

  /*
   * Here we check for PyDatum, a wrapper class around Datum, defined in the cython module kernel.pyx.
   * and made available in "kernel.h"
   */ 
  if (PyObject_TypeCheck(pObj, &PyTokenType))
  { // Object is encapsulated Datum
      Token* t = reinterpret_cast<PyToken*>(pObj)->thisptr;
      if( not t)
      {
          const char error[]="PyToken must contain a valid Token.\n";
          PyErr_SetString(NESTError_, error);
          return 0;
      }

      Datum *d=t->datum();
      if(d ) // In case pObj was uninitialized
      {
          d->addReference();
          return d;
      }
      else
      {
          const char error[]="PyToken must contain a valid Datum reference.\n";
          PyErr_SetString(NESTError_, error);
          return 0;
      }
  }

#if PY_MAJOR_VERSION >= 3
  if (PyUnicode_Check(pObj)) // object is string
  {
      PyObject *byte_repr=PyUnicode_AsUTF8String(pObj); // We re-code the unicode string into a bytes object
      if(byte_repr ==0)
      {
	  PyErr_SetString(NESTError_, "could not encode unicode as 'UTF-8'.");
	  return 0;
      }
      StringDatum *result=new StringDatum(PyBytes_AsString(byte_repr));
      Py_DECREF(byte_repr);
      return result;
  }
#endif

#ifdef HAVE_NUMPY
  if (PyArray_CheckScalar(pObj)) { // handle numpy array scalars

    PyArray_Descr *typecode;
    typecode = PyArray_DescrFromScalar(pObj);

    switch (typecode->type_num)
    {
      case NPY_INT:
      {
        int val;
        PyArray_ScalarAsCtype(pObj, &val);
        return new IntegerDatum(val);
      }
      case NPY_LONG:
      {
        long val;
        PyArray_ScalarAsCtype(pObj, &val);
        return new IntegerDatum(val);
      }
      case NPY_DOUBLE:
      {
        double val;
        PyArray_ScalarAsCtype(pObj, &val);
        return new DoubleDatum(val);
      }
      case NPY_OBJECT: // handle 0-dim numpy array, which are treated as scalars
      {
        PyArrayObject *array = 0;
        array = (PyArrayObject*) pObj;
        assert(array != 0);
        return PyObject_as_Datum(PyArray_ToScalar(array->data, pObj));
      }
      default:
      {
        std::string error = String::compose("Unsupported Numpy array scalar type: '%1'.\n"
                                            "If you think this is an error, tell us at nest_user@nest-initiative.org",
                                            typecode->type_num);
        PyErr_SetString(NESTError_, error.c_str());
        return NULL;
      }
    }
  }

  if (PyArray_Check(pObj)) { // handle numpy arrays by sending as VectorDatum
    
    PyArrayObject *array = 0;
    array = (PyArrayObject*) pObj;
    assert(array != 0);

    Datum* vd = 0;
    int size = PyArray_Size(pObj);

    // raise an exception if array's dimensionality is not 1
    if (array->nd != 1)
    {
      std::string error = String::compose("The given Numpy array has an unsupported dimensionality of %1.\n"
                                          "Only one-dimensional arrays are supported. "
                                          "If you think this is an error,\ntell us at nest_user@nest-initiative.org",
                                          array->nd);
      PyErr_SetString(NESTError_, error.c_str());
      return NULL;
    }

    switch (array->descr->type_num)
    {
      case NPY_INT :
      case NPY_LONG:
      {
        long *begin = reinterpret_cast<long*>(array->data);
        std::vector<long>* datavec;

        // Check if we really have a pure 1-dim array instead of a selection like obtained with a[:,1]
        // The latter needs to be treated different because we have to move in steps of array->strides[0]
        if (array->strides[0] == sizeof(long))
          datavec = new std::vector<long>(begin, begin + size);
        else
        {
          datavec = new std::vector<long>(size);
          for (int i = 0; i < size; i++)
            (*datavec)[i] = *(long*)(array->data + i*array->strides[0]);
        }

        vd = new IntVectorDatum(datavec);
        break;
      }
      case NPY_DOUBLE:
      {
        double *begin = reinterpret_cast<double*>(array->data);
        std::vector<double>* datavec;
        
        // Check if we really have a pure 1-dim array instead of a selection like obtained with a[:,1]
        // The latter needs to be treated different because we have to move in steps of array->strides[0]
        if (array->strides[0] == sizeof(double))
          datavec = new std::vector<double>(begin, begin + size);
        else
        {
          datavec = new std::vector<double>(size);
          for (int i = 0; i < size; i++)
            (*datavec)[i] = *(double*)(array->data + i*array->strides[0]);
        }

        vd = new DoubleVectorDatum(datavec);
        break;
      }
      default:
      {
        std::string error = String::compose("Unsupported Numpy array type: '%1'.\n"
                                            "If you think this is an error, tell us at bugs@nest-initiative.org",
                                            array->descr->type_num);
        PyErr_SetString(NESTError_, error.c_str());
        return NULL;
      }
    }
    
    return vd;
  }
#endif

  if (PyList_Check(pObj) || PyTuple_Check(pObj)) { // object is a list or a tuple
    ArrayDatum *d = new ArrayDatum();
    PyObject* subPyObj;

    size_t size = (PyList_Check(pObj)) ? PyList_Size(pObj) : PyTuple_Size(pObj);

    d->reserve(size);
 
    for (size_t i = 0; i < size; ++i) {
      subPyObj = (PyList_Check(pObj)) ? PyList_GetItem(pObj, i) : PyTuple_GetItem(pObj, i);
      Datum *tmp_d=PyObject_as_Datum(subPyObj);
      assert(tmp_d !=0);
      d->push_back(tmp_d);
      if (PyErr_Occurred()) {
        delete d;
        return NULL;
      }
    }
    return d;
  }

  if (PyDict_Check(pObj)) { // object is a dictionary
    DictionaryDatum *d = new DictionaryDatum(new Dictionary());

    PyObject* subPyObj=0;
    PyObject *key=0;
    Py_ssize_t pos = 0;

    while (PyDict_Next(pObj, &pos, &key, &subPyObj)) 
    {
      Token t(PyObject_as_Datum(subPyObj));
      if (PyErr_Occurred()) {
	delete d;
	return NULL;
      }
	 
#if PY_MAJOR_VERSION < 3
     if (PyString_Check(key)) 
	{
	  (*d)->insert(PyString_AsString(key), t);
	}
#else
	if (PyUnicode_Check(key)) 
	{
	    PyObject *byte_repr=PyUnicode_AsUTF8String(key); // We re-code the unicode string into a bytes object
	    if(byte_repr ==0)
	    {
		PyErr_SetString(NESTError_, "could not encode unicode as 'UTF-8' Skipping key.");
		continue;
	    }
	    (*d)->insert(PyBytes_AsString(byte_repr), t);
	    Py_DECREF(byte_repr);
	}
#endif
	else
	{ 
	    PyObject* keystr = PyObject_Str(key);
	    if (keystr == NULL)
	    {
		PyErr_Warn(PyExc_Warning, "Non-string key in Python dictionary. Using bogus key in  dictionary.");
		(*d)->insert("BOGUS_KEY", t);
	    }
	    else
		(*d)->insert(PyString_AsString(pObj), t);
	}
    }
    return d;
  }

#if PY_MAJOR_VERSION >= 3
  if(PySequence_Check(pObj)) { // object is a sequence

     PyObject* list = PySequence_List(pObj);

     Py_DECREF(pObj);

    return PyObject_as_Datum(list);
  }
# endif

  return new PyObjectDatum(pObj);


  std::string error = String::compose("Python object of type '%1' cannot be converted to SLI.\n"
                                      "If you think this is an error, tell us at nest_user@nest-initiative.org",
                                      pObj->ob_type->tp_name);
  PyErr_SetString(NESTError_, error.c_str());
  return 0;
}


