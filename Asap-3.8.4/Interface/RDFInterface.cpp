// -*- C++ -*-
// RDFInterface.cpp: Python interface to the RawRadialDistribution function.
//
// Copyright (C) 2008 Jakob Schiotz and Center for Individual
// Nanoparticle Functionality, Department of Physics, Technical
// University of Denmark.  Email: schiotz@fysik.dtu.dk
//
// This file is part of Asap version 3.
//
// This program is free software: you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License
// version 3 as published by the Free Software Foundation.  Permission
// to use other versions of the GNU Lesser General Public License may
// granted by Jakob Schiotz or the head of department of the
// Department of Physics, Technical University of Denmark, as
// described in section 14 of the GNU General Public License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// and the GNU Lesser Public License along with this program.  If not,
// see <http://www.gnu.org/licenses/>.

#include "RDFInterface.h"
#include "RawRadialDistribution.h"
#include "ExceptionInterface.h"
#include "PythonConversions.h"
#include <set>
using std::set;
#include <iostream>

namespace ASAPSPACE {
PyObject *PyAsap_RawRDF(PyObject *noself, PyObject *args)
{
  // Arguments are atoms, rMax, nBins, groups, ngroups, elements

  PyObject *atoms;
  double rMax;
  long nBins_arg;
  PyArrayObject *groups;
  int nGroups;
  PyObject *pyelements;
  if (!PyArg_ParseTuple(args, "OdlO!iO:RawRDF", &atoms, &rMax, &nBins_arg,
			&PyArray_Type, &groups, &nGroups, &pyelements))
    return NULL;
  npy_intp nBins = nBins_arg;
  // Check that the groups array is sensible.
  if (PyArray_NDIM(groups) != 1           // one-dimensional
      || PyArray_TYPE(groups) != NPY_INT32   // array of int
      || !PyArray_ISCARRAY_RO(groups))
    {
      PyErr_SetString(PyExc_ValueError,
	  "RawRDF: groups must be a contiguous array of nAtoms Int32s");
      return NULL;
    }
  
  // Make parallel C++ and Python data structures for the RDFs,
  // sharing memory so the Python objects get filled when the C++
  // objects are updated.

  // "result" is the Python tuple containing the result values from
  // this function.  It owns all subsequently created Python
  // objects, so if an error occurs, deallocating this one should
  // release all memory.
  PyObject *result = PyTuple_New(3);
  if (!result) 
    return NULL;

  // Create the data structures, being prepared to release 'result'
  // if an error occurs.
  int errln = 0;  // Where did the error occur
#define CHK(x) if (!x) {errln = __LINE__; throw AsapError("Oops");}
#define CHKINT(x) if (x<0) {errln = __LINE__; throw AsapError("Oops");}
  PyObject *globalRdf = NULL;
  PyObject *rdfPy = NULL;
  PyObject *rdfcountPy = NULL;
  RDFtype rdf(nGroups);
  RDFcountType rdfcount(nGroups);
  try
    {
      // The global RDF array
      globalRdf = PyArray_ZEROS(1, &nBins, NPY_LONG, 0);
      CHK(globalRdf);
      CHKINT(PyTuple_SetItem(result, 0, globalRdf));
      // The main RDF data structure
      rdfPy = PyList_New(nGroups);
      CHK(rdfPy);
      CHKINT(PyTuple_SetItem(result, 1, rdfPy));
      // The structure for counting atoms in each sub-RDF
      rdfcountPy = PyList_New(nGroups);
      CHK(rdfcountPy);
      CHKINT(PyTuple_SetItem(result, 2, rdfcountPy));
      // Now populate both Python lists
      for (int i = 0; i < nGroups; i++)
	{
	  PyObject *tmp;
	  tmp = PyDict_New();
	  CHK(tmp);
	  CHKINT(PyList_SET_ITEM(rdfPy, i, tmp));
	  tmp = PyDict_New();
	  CHK(tmp);
	  CHKINT(PyList_SET_ITEM(rdfcountPy, i, tmp));
	}
      // Now we have two lists of empty dictionaries in each
      // language.  The Python RDF dictionary should now be
      // populated with pairs of elements mapping to NumPy
      // arrays containing the RDF.  The C++ maps should be
      // populated with the same pairs mapping to pointers into
      // these arrays.
      set<int> elements;
      CHKINT(PyAsap_SetIntFromArray(elements, pyelements));
      for (int i = 0; i < nGroups; i++)
	{
	  PyObject *mapPy = PyList_GetItem(rdfPy, i);
	  CHK(mapPy);
	  for (set<int>::const_iterator e1 = elements.begin();
	       e1 != elements.end(); ++e1)
	    for (set<int>::const_iterator e2 = elements.begin();
		 e2 != elements.end(); ++e2)
	      {
		// Fill the Python dict
		PyObject *keyPy = NULL;
		PyObject *valPy = NULL;
		try
		  {
		    keyPy = Py_BuildValue("(ii)", *e1, *e2);
		    CHK(keyPy);

		    valPy = PyArray_ZEROS(1, &nBins, NPY_LONG, 0);
		    CHK(valPy);

		    CHKINT(PyDict_SetItem(mapPy, keyPy, valPy));
		    Py_DECREF(keyPy); // mapPy now owns it.
		    Py_DECREF(valPy);
		  }
		catch(AsapError ex)
		  {
		    Py_XDECREF(keyPy);
		    Py_XDECREF(valPy);
		    throw;
		  }
		// Fill the C++ map
		pair<int, int> key(*e1, *e2);
		rdf[i][key] = (long *) PyArray_DATA((PyArrayObject *) valPy);
	      }
	}
    }
  catch(AsapError ex)
    {
      Py_DECREF(result);  // Discard all objects created so far
      std::cerr << "Error in RadialDistributionFunction (memory?) ("
	   << __FILE__ << ":" << errln << ")" << std::endl;
      if (!PyErr_Occurred())
	PyErr_Format(PyExc_MemoryError,
		     "Out of memory in RadialDistributionFunction (%s:%d)",
		     __FILE__, errln);
      return NULL;
    }
  // Now the data structures are ready.  Call the function doing
  // the actual work, but be prepared to deallocate the return
  // value if an error occurs.
  try
    {
      CHECKNOASAPERROR;
      RawRadialDistribution(atoms, nGroups, (int *) PyArray_DATA(groups),
			    rMax, (int) nBins, rdf, rdfcount,
			    (long *) PyArray_DATA((PyArrayObject *) globalRdf));
      PROPAGATEASAPERROR;
    }
  catch(AsapError ex)
    {
      Py_DECREF(result);
      string msg = ex.GetMessage();
      PyErr_SetString(PyAsap_ErrorObject, msg.c_str());
      return NULL;
    }
  // The result has automagically been inserted into the Python
  // objects except for rdfcount.  Insert that now.
  try
    {
      for (int i = 0; i < nGroups; i++)
	{
	  PyObject *mapPy = PyList_GetItem(rdfcountPy, i);
	  CHK(mapPy);
	  assert(PyDict_Size(mapPy) == 0);
	  map<int,long>::const_iterator enditer = rdfcount[i].end();
	  for (map<int,long>::const_iterator mapiter = rdfcount[i].begin();
	       mapiter != enditer; ++mapiter)
	    {
	      PyObject *key = PyInt_FromLong(mapiter->first);
	      CHK(key);
	      PyObject *val = PyInt_FromLong(mapiter->second);
	      CHK(val);
	      int x = PyDict_SetItem(mapPy, key, val);
	      Py_DECREF(key);
	      Py_DECREF(val);
	      CHKINT(x);
	    }
	}
    } 
  catch(AsapError ex)
    {
      Py_DECREF(result);  // Discard all objects created so far
      std::cerr << "Error in RadialDistributionFunction (memory?) ("
	   << __FILE__ << ":" << errln << ")" << std::endl;
      if (!PyErr_Occurred())
	PyErr_Format(PyExc_RuntimeError,
		     "Unknown error in RadialDistributionFunction (%s:%d)",
		     __FILE__, errln);
      return NULL;
    }
  return result;
  // Undefine the macros so we get an error in the unlikely case the
  // name was already in use.
}

} // end namespace

