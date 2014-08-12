#ifndef PY3K_COMPATIBILITY_H
#define  PY3K_COMPATIBILITY_H

#if PY_MAJOR_VERSION >= 3
#define IS_PY3K
#define PyInt_Check PyLong_Check
#define PyInt_AsLong PyLong_AsLong
#define PyInt_FromLong PyLong_FromLong
#define PyString_Check PyBytes_Check
#define PyString_AsString PyBytes_AsString
#define PyString_FromString PyBytes_FromString
#endif

#ifndef DL_IMPORT
#define DL_IMPORT(t) t
#endif

#endif
