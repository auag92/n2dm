// -*- C++ -*-
// PotentialInterface.cpp: Python interface to the calculator objects.
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

#define MCEMT

#include "AsapPython.h"
#include "Asap.h"
#include "PotentialInterface.h"
#include "ExceptionInterface.h"
#include "PythonConversions.h"
#include "Templates.h"
#include "EMT.h"
#include "EMT2013.h"
#include "MonteCarloEMT.h"
#include "RGL.h"
#include "LennardJones.h"
#include "Morse.h"
#include "BrennerPotential.h"
#include "EMTParameterProviderInterface.h"
#include "ImagePotential.h"
//#define ASAPDEBUG
#include "Debug.h"

// PyAsap_PotentialObject is defined in Potential.h

// Potential base class
namespace ASAPSPACE {
PyTypeObject PyAsap_PotentialType = {
  PyObject_HEAD_INIT(NULL)
  0,
  "asapserial3.Potential",
  sizeof(PyAsap_PotentialObject),
  // The rest are initialized by name for reliability.
};
} // namespace

static char Potential_Docstring[] = "Asap Potential abstract base class.\n";


static PyTypeObject PyAsap_EMTType = {
  PyObject_HEAD_INIT(NULL)
  0,
  "asapserial3.EMT",
  sizeof(PyAsap_PotentialObject),
  // The rest are initialized by name for reliability.
};
  
static char EMT_Docstring[] = "Effective Medium Theory calculator.\n";

static PyTypeObject PyAsap_EMT2013Type = {
  PyObject_HEAD_INIT(NULL)
  0,
  "asapserial3.EMT2013",
  sizeof(PyAsap_PotentialObject),
  // The rest are initialized by name for reliability.
};

static char EMT2013_Docstring[] = "Effective Medium Theory version 2011 calculator.\n";

#ifdef MCEMT
static PyTypeObject PyAsap_MonteCarloEMTType = {
  PyObject_HEAD_INIT(NULL)
  0,
  "asapserial3.MonteCarloEMT",
  sizeof(PyAsap_PotentialObject),
  // The rest are initialized by name for reliability.
};
  
static char MonteCarloEMT_Docstring[] = "Effective Medium Theory calculator optimized for Monte Carlo simulations.\n";
#endif

static PyTypeObject PyAsap_RGLType = {
  PyObject_HEAD_INIT(NULL)
  0,
  "asapserial3.RGL",
  sizeof(PyAsap_PotentialObject),
  // The rest are initialized by name for reliability.
};
  
static char RGL_Docstring[] = "RGL tight-binding potential calculator.\n";

static PyTypeObject PyAsap_LennardJonesType = {
  PyObject_HEAD_INIT(NULL)
  0,
  "asapserial3.LennardJones",
  sizeof(PyAsap_PotentialObject),
  // The rest are initialized by name for reliability.
};
  
static char LennardJones_Docstring[] = "Lennard-Jones calculator.\n";

static PyTypeObject PyAsap_MorseType = {
  PyObject_HEAD_INIT(NULL)
  0,
  "asapserial3.Morse",
  sizeof(PyAsap_PotentialObject),
  // The rest are initialized by name for reliability.
};
  
static char Morse_Docstring[] = "Morse potential calculator.\n";

static PyTypeObject PyAsap_BrennerType = {
  PyObject_HEAD_INIT(NULL)
  0,
  "asapserial3.BrennerPotential",
  sizeof(PyAsap_PotentialObject),
  // The rest are initialized by name for reliability.
};
  
static char Brenner_Docstring[] =
  "Brenner potential calculator for C, H, Si and Ge.\n";

// A few convenience macros

#define CHECK_POT_UNINIT if (self->cobj != NULL) {		\
    PyErr_SetString(PyAsap_ErrorObject,				\
		    "Potential object already initialized.");	\
    return -1;							\
  }

#define CHECK_POT_INIT if (self->cobj == NULL) {	  \
    PyErr_SetString(PyAsap_ErrorObject,			  \
		    "Potential object not initialized."); \
    return NULL;					  \
  }

namespace ASAPSPACE {
int PyAsap_PotentialInit(PyAsap_PotentialObject *self, PyObject *args,
                         PyObject *kwargs)
{
  CHECK_POT_UNINIT;
  self->weakrefs = NULL;
  self->setatoms_called = false;
  return 0;
}
} // namespace

static int PyAsap_EMTInit(PyAsap_PotentialObject *self, PyObject *args,
			  PyObject *kwargs)
{
  static char *kwlist[] = {"prov", NULL};
  
  PyObject *provider = NULL;
  if (!PyArg_ParseTupleAndKeywords(args, kwargs,  "O!", kwlist,
				   &PyAsap_EMTParamProvType, &provider))
    return -1;
  assert(provider != NULL);
  if (PyAsap_PotentialType.tp_init((PyObject *)self, args, kwargs) < 0)
    return -1;
  try
    {
      self->cobj = new EMT((PyObject *)self, provider);
      self->orig_cobj = self->cobj;
    }
  catch (AsapError &e)
    {
      string msg = e.GetMessage();
      PyErr_SetString(PyAsap_ErrorObject, msg.c_str());
      return -1;
    }
  catch (AsapPythonError &e)
    {
      return -1;
    }
  if (self->cobj == NULL)
    return -1;
  return 0;
}

static int PyAsap_EMT2013Init(PyAsap_PotentialObject *self, PyObject *args,
                              PyObject *kwargs)
{
  static char *kwlist[] = {"parameters", "no_new_elements", NULL};

  PyObject *parameters = NULL;
  int no_new = 0;
  if (!PyArg_ParseTupleAndKeywords(args, kwargs,  "O!i", kwlist,
                                   &PyDict_Type, &parameters, &no_new))
    return -1;
  assert(parameters != NULL);
  bool no_new_parameters = (bool) no_new;
  if (PyAsap_PotentialType.tp_init((PyObject *)self, args, kwargs) < 0)
    return -1;
  try
    {
      self->cobj = new EMT2013((PyObject *)self, parameters, no_new_parameters);
      self->orig_cobj = self->cobj;
    }
  catch (AsapError &e)
    {
      string msg = e.GetMessage();
      PyErr_SetString(PyAsap_ErrorObject, msg.c_str());
      return -1;
    }
  catch (AsapPythonError &e)
    {
      return -1;
    }
  if (self->cobj == NULL)
    return -1;
  return 0;
}

#ifdef MCEMT
static int PyAsap_MonteCarloEMTInit(PyAsap_PotentialObject *self,
				    PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"prov", NULL};
  
  PyObject *provider = NULL;
  if (!PyArg_ParseTupleAndKeywords(args, kwargs,  "O!", kwlist,
				   &PyAsap_EMTParamProvType, &provider))
    return -1;
  assert(provider != NULL);
  if (PyAsap_PotentialType.tp_init((PyObject *)self, args, kwargs) < 0)
    return -1;
  self->cobj = new MonteCarloEMT((PyObject *)self, provider);
  self->orig_cobj = self->cobj;
  if (self->cobj == NULL)
    return -1;
  return 0;
}
#endif

static int PyAsap_RGLInit(PyAsap_PotentialObject *self, PyObject *args,
			    PyObject *kwargs)
{
  static char *kwlist[] = {"elements", "p", "q", "A", "qsi2", "r0",
                           "p3", "p4", "p5", "q3", "q4", "q5",
                           "rcs", "rce", NULL};
  
  PyObject *elements_obj;
  PyObject *p_obj;
  PyObject *q_obj;
  PyObject *A_obj;
  PyObject *qsi2_obj;
  PyObject *r0_obj;
  PyObject *p3_obj;
  PyObject *p4_obj;
  PyObject *p5_obj;
  PyObject *q3_obj;
  PyObject *q4_obj;
  PyObject *q5_obj;
  double rcs;
  double rce;
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOOOOOOOOOOdd", kwlist,
				                           &elements_obj, &p_obj, &q_obj, &A_obj,
                                   &qsi2_obj, &r0_obj, &p3_obj, &p4_obj, &p5_obj,
                                   &q3_obj, &q4_obj, &q5_obj, &rcs, &rce))
    return -1;
  self->weakrefs = NULL;
  if (PyAsap_PotentialType.tp_init((PyObject *)self, args, kwargs) < 0)
    return -1;
  vector<int> elements;
  TinyMatrix<double> p;
  TinyMatrix<double> q;
  TinyMatrix<double> A;
  TinyMatrix<double> qsi2;
  TinyMatrix<double> r0;
  TinyMatrix<double> p3;
  TinyMatrix<double> p4;
  TinyMatrix<double> p5;
  TinyMatrix<double> q3;
  TinyMatrix<double> q4;
  TinyMatrix<double> q5;
  if (PyAsap_VectorIntFromArray(elements, elements_obj) ||
      PyAsap_TinyMatrixDoubleFromArray(p, p_obj) ||
      PyAsap_TinyMatrixDoubleFromArray(q, q_obj) ||
      PyAsap_TinyMatrixDoubleFromArray(A, A_obj) ||
      PyAsap_TinyMatrixDoubleFromArray(qsi2, qsi2_obj) ||
      PyAsap_TinyMatrixDoubleFromArray(r0, r0_obj) ||
      PyAsap_TinyMatrixDoubleFromArray(p3, p3_obj) ||
      PyAsap_TinyMatrixDoubleFromArray(p4, p4_obj) ||
      PyAsap_TinyMatrixDoubleFromArray(p5, p5_obj) ||
      PyAsap_TinyMatrixDoubleFromArray(q3, q3_obj) ||
      PyAsap_TinyMatrixDoubleFromArray(q4, q4_obj) ||
      PyAsap_TinyMatrixDoubleFromArray(q5, q5_obj))
    return -1;
  self->cobj = new RGL((PyObject *)self, elements, p, q, A, qsi2, r0, p3, p4,
                       p5, q3, q4, q5, rcs, rce);
  self->orig_cobj = self->cobj;
  if (self->cobj == NULL)
    return -1;
  return 0;
}

static int PyAsap_LennardJonesInit(PyAsap_PotentialObject *self, PyObject *args,
				   PyObject *kwargs)
{
  static char *kwlist[] = {"numElements", "elements", "epsilon", "sigma",
			   "masses", "rCut", "modified", NULL};
  
  int numElements;
  PyObject *elements_obj;
  PyObject *epsilon_obj;
  PyObject *sigma_obj;
  PyObject *masses_obj;
  double rCut = -1.0;
  char modified = 1;
  if (!PyArg_ParseTupleAndKeywords(args, kwargs,  "iOOOO|db", kwlist,
				   &numElements, &elements_obj, &epsilon_obj,
				   &sigma_obj, &masses_obj, &rCut, &modified))
    return -1;
  if (PyAsap_PotentialType.tp_init((PyObject *)self, args, kwargs) < 0)
    return -1;
  vector<int> elements;
  vector<double> epsilon;
  vector<double> sigma;
  vector<double> masses;
  if (PyAsap_VectorIntFromArray(elements, elements_obj)
      || PyAsap_VectorDoubleFromArray(epsilon, epsilon_obj)
      || PyAsap_VectorDoubleFromArray(sigma, sigma_obj)
      || PyAsap_VectorDoubleFromArray(masses, masses_obj))
    return -1;
  self->cobj = new LennardJones((PyObject *)self,
                                numElements, elements, epsilon, sigma, masses,
			        rCut, (bool)modified);
  self->orig_cobj = self->cobj;
  if (self->cobj == NULL)
    return -1;
  return 0;
}

static int PyAsap_MorseInit(PyAsap_PotentialObject *self, PyObject *args,
			    PyObject *kwargs)
{
  static char *kwlist[] = {"elements", "epsilon", "alpha", "rmin",
			   "rCut", "modified", NULL};
  
  PyObject *elements_obj;
  PyObject *epsilon_obj;
  PyObject *alpha_obj;
  PyObject *rmin_obj;
  double rCut = -1.0;
  char modified = 1;
  if (!PyArg_ParseTupleAndKeywords(args, kwargs,  "OOOO|db", kwlist,
				   &elements_obj, &epsilon_obj,
				   &alpha_obj, &rmin_obj, &rCut, &modified))
    return -1;
  if (PyAsap_PotentialType.tp_init((PyObject *)self, args, kwargs) < 0)
    return -1;
  vector<int> elements;
  vector<double> epsilon;
  vector<double> alpha;
  vector<double> rmin;
  if (PyAsap_VectorIntFromArray(elements, elements_obj)
      || PyAsap_VectorDoubleFromArray(epsilon, epsilon_obj)
      || PyAsap_VectorDoubleFromArray(alpha, alpha_obj)
      || PyAsap_VectorDoubleFromArray(rmin, rmin_obj))
    return -1;
  self->cobj = new Morse((PyObject *)self, elements, epsilon, alpha, rmin,
			 rCut, (bool)modified);
  self->orig_cobj = self->cobj;
  if (self->cobj == NULL)
    return -1;
  return 0;
}

static int PyAsap_BrennerInit(PyAsap_PotentialObject *self,
			      PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {NULL};
  
  if (!PyArg_ParseTupleAndKeywords(args, kwargs,  "", kwlist))
    return -1;
  if (PyAsap_PotentialType.tp_init((PyObject *)self, args, kwargs) < 0)
    return -1;
  self->cobj = new BrennerPotential((PyObject *)self);
  self->orig_cobj = self->cobj;
  if (self->cobj == NULL)
    return -1;
  return 0;
}

static void PyAsap_PotentialDealloc(PyObject *self)
{
  PyAsap_PotentialObject *myself = (PyAsap_PotentialObject *) self;
  if ((myself->cobj != NULL) && (myself->cobj != myself->orig_cobj))
    delete myself->cobj;       // Delete wrapper if separate object.
  if (myself->orig_cobj != NULL)
    delete myself->orig_cobj;  // Delete the actual Potential object.
  self->ob_type->tp_free(self);
}


static PyObject *PyAsap_PotentialGetPotentialEnergy(PyAsap_PotentialObject *self,
					      PyObject *args)
{
  PyObject *atoms = NULL;
  if (!PyArg_ParseTuple(args, "O", &atoms))
    return NULL;
  CHECK_POT_INIT;
  double e;
  try {
    e = self->cobj->GetPotentialEnergy(atoms);
  }
  POTCATCHEXCEPTION;
  return PyFloat_FromDouble(e);
}

static PyObject *PyAsap_PotentialSetAtoms(PyAsap_PotentialObject *self, PyObject *args)
{
  PyObject *atoms = NULL;
  PyObject *py_accessobj = NULL;
  if (!PyArg_ParseTuple(args, "O|O", &atoms, &py_accessobj))
    return NULL;
  CHECK_POT_INIT;
  Atoms *accessobj = NULL;
  if ((py_accessobj != NULL) && (py_accessobj != Py_None))
    {
      accessobj = (Atoms *) PyCObject_AsVoidPtr(py_accessobj);
    }
  try {
#if 0
      std::cerr << "Potential.SetAtoms: cobj: "
          << self->cobj << " " << self->cobj->GetName()
          << "     orig_cobj: "<< self->orig_cobj << " " << self->orig_cobj->GetName() << std::endl;
#endif
      self->cobj->SetAtoms(atoms, accessobj);
  }
  POTCATCHEXCEPTION;
  self->setatoms_called = true;
  Py_RETURN_NONE;
}

static PyObject *PyAsap_PotentialGetPotentialEnergies(PyAsap_PotentialObject *self,
						PyObject *args)
{
  PyObject *atoms = NULL;
  if (!PyArg_ParseTuple(args, "O", &atoms))
    return NULL;
  CHECK_POT_INIT;
  try {
    const vector<double> &energies = self->cobj->GetPotentialEnergies(atoms);
    return PyAsap_ArrayFromVectorDouble(energies);
  }
  POTCATCHEXCEPTION;
}

static PyObject *PyAsap_PotentialGetForces(PyAsap_PotentialObject *self,
				     PyObject *args)
{
  PyObject *atoms = NULL;
  if (!PyArg_ParseTuple(args, "O", &atoms))
    return NULL;
  CHECK_POT_INIT;
  try {
    FP_EXCEPT_ON;
    const vector<Vec> &forces = self->cobj->GetForces(atoms);
    FP_EXCEPT_OFF;
    return PyAsap_ArrayFromVectorVec(forces);
  }
  POTCATCHEXCEPTION;
}

static PyObject *PyAsap_PotentialGetStress(PyAsap_PotentialObject *self,
				     PyObject *args)
{
  PyObject *atoms = NULL;
  if (!PyArg_ParseTuple(args, "O", &atoms))
    return NULL;
  CHECK_POT_INIT;
  vector<double> stress(6);
  try {
      SymTensor s = self->cobj->GetStress(atoms);
      for (int i = 0; i < 6; i++)
        stress[i] = s[i];  // Should instead write a PyAsap_ArrayFromSymTensor
  }
  POTCATCHEXCEPTION;
  return PyAsap_ArrayFromVectorDouble(stress);
}

static PyObject *PyAsap_PotentialGetStresses(PyAsap_PotentialObject *self,
					     PyObject *args)
{
  PyObject *atoms = NULL;
  if (!PyArg_ParseTuple(args, "O", &atoms))
    return NULL;
  CHECK_POT_INIT;
  try {
    const vector<SymTensor> &stresses = self->cobj->GetStresses(atoms);
    return PyAsap_ArrayFromVectorSymTensor(stresses);
  }
  POTCATCHEXCEPTION;
}

static PyObject *PyAsap_PotentialGetCutoff(PyAsap_PotentialObject *self,
                                           PyObject *noargs)
{
  CHECK_POT_INIT;
  double cutoff;
  try {
    cutoff = self->cobj->GetCutoffRadius();
  }
  POTCATCHEXCEPTION;
  return PyFloat_FromDouble(cutoff);
}

static PyObject *PyAsap_PotentialGetExtra(PyAsap_PotentialObject *self,
					  PyObject *args)
{
  const char *property;
  if (!PyArg_ParseTuple(args, "s", &property))
    return NULL;
  CHECK_POT_INIT;
  EMT2013 *emt2013 = dynamic_cast<EMT2013 *>(self->orig_cobj);
  if (emt2013 != NULL && strcmp(property, "parameters") == 0)
    {
      return emt2013->GetParameterDict();
    }
  EMT *emt = dynamic_cast<EMT*>(self->orig_cobj);
  if (emt != NULL && strcmp(property, "sigma") == 0)
    {
      const vector<vector<double> > &s1 = emt->GetSigma1();
      const vector<vector<double> > &s2 = emt->GetSigma2();
      PyObject *py_s1 = PyList_New(s1.size());
      PyObject *py_s2 = PyList_New(s2.size());
      assert(s1.size() == s2.size());
      for (int i = 0; i < s1.size(); i++)
	{
	  PyList_SET_ITEM(py_s1, i, PyAsap_ArrayFromVectorDouble(s1[i]));
	  PyList_SET_ITEM(py_s2, i, PyAsap_ArrayFromVectorDouble(s2[i]));
	}
      return Py_BuildValue("NN", py_s1, py_s2);
    }
  PyErr_SetString(PyExc_ValueError, "Unknown extra property for this potential");
  return NULL;
}

static PyObject *PyAsap_EMTSetSubtractE0(PyAsap_PotentialObject *self,
                                          PyObject *args)
{
  int subtractE0;
  if (!PyArg_ParseTuple(args, "i", &subtractE0))
    return NULL;
  CHECK_POT_INIT;
  EMT *emt = dynamic_cast<EMT*>(self->orig_cobj);
  assert(emt != NULL);  // Should not be possible
  emt->SetSubtractE0(subtractE0);
  Py_RETURN_NONE;
}

// Four static strings used when comparing.
static PyObject *PyAsap_String_Energy;
static PyObject *PyAsap_String_Forces;
static PyObject *PyAsap_String_Stress;

static PyObject *PyAsap_CalcReq(PyAsap_PotentialObject *self,
				PyObject *args)
{
  PyObject *atoms = NULL;
  PyObject *proplist = NULL;
  if (!PyArg_ParseTuple(args, "OO:calculation_required", &atoms, &proplist))
    return NULL;
  assert(proplist != NULL);
  if (!PySequence_Check(proplist))
    {
      PyErr_SetString(PyExc_TypeError, "Argument to calculation_required must be a sequence");
      return NULL;
    }
  bool result = false;
  Py_ssize_t n = PySequence_Size(proplist); 
  for (Py_ssize_t i = 0; i < n; i++)
    {
      PyObject *pystr = PySequence_GetItem(proplist, i);
      assert(pystr != NULL);
      if (!PyString_Check(pystr))
	{
	  Py_DECREF(pystr);
	  PyErr_SetString(PyExc_ValueError, "Non-string passed to calculation_required.");
	  return NULL;
	}
      if (PyObject_Compare(pystr, PyAsap_String_Energy) == 0)
	result = result || self->cobj->CalcReq_Energy(atoms);
      else if (PyObject_Compare(pystr, PyAsap_String_Forces) == 0)
	result = result || self->cobj->CalcReq_Forces(atoms);
      else if (PyObject_Compare(pystr, PyAsap_String_Stress) == 0)
	result = result || self->cobj->CalcReq_Stress(atoms);
      else 
        result = true;  // All unsupported keywords must cause True to be returned
      Py_DECREF(pystr);
    }
  if (result)
    Py_RETURN_TRUE;
  else
    Py_RETURN_FALSE;
}

static PyObject *PyAsap_PotentialGetNbList(PyAsap_PotentialObject *self,
					   PyObject *noargs)
{
  CHECK_POT_INIT;
  PyObject *nbl = self->cobj->GetNeighborList();
  if (nbl == NULL)
    PyErr_SetString(PyExc_RuntimeError,
		    "No neighbor list (potential still unused?)");
  else
    Py_INCREF(nbl);
  return nbl;
}

namespace ASAPSPACE {

char PyAsap_PotentialUseImageAtoms_Docstring[] =
    "Use image atoms (disable minimum image convention). INTERNAL USE ONLY!";

PyObject *PyAsap_PotentialUseImageAtoms(PyAsap_PotentialObject *self,
                                        PyObject *noargs)
{
  CHECK_POT_INIT;
  if (self->cobj != self->orig_cobj)
    return PyErr_Format(PyExc_RuntimeError,
        "Error: _use_imageatoms called, but %s object is already wrapped by a %s object.",
        self->orig_cobj->GetName().c_str(),
        self->cobj->GetName().c_str() );
  if (self->setatoms_called)
    return PyErr_Format(PyExc_RuntimeError,
        "Error: _use_imageatoms called, but %s object has already seen the atoms.",
        self->orig_cobj->GetName().c_str());
  try
    {
      DEBUGPRINT;
      self->cobj = new ImagePotential((PyObject *)self, self->orig_cobj);
      DEBUGPRINT;
    }
  catch (AsapError &e)
    {
      string msg = e.GetMessage();
      PyErr_SetString(PyAsap_ErrorObject, msg.c_str());
      return NULL;
    }
  catch (AsapPythonError &e)
    {
      return NULL;
    }
  Py_RETURN_NONE;
}

} // namespace

static PyObject *PyAsap_Parallel(PyAsap_PotentialObject *self,
				 PyObject *noargs)
{
  CHECK_POT_INIT;
  if (self->cobj->Parallelizable())
    Py_RETURN_TRUE;
  else
    Py_RETURN_FALSE;
}

static PyObject *PyAsap_PotPrintMemory(PyAsap_PotentialObject *self,
				    PyObject *noargs)
{
  CHECK_POT_INIT;
  long mem = self->cobj->PrintMemory();
  return PyInt_FromLong(mem);
}

static PyMethodDef PyAsap_PotentialMethods[] = {
  {"set_atoms", (PyCFunction)PyAsap_PotentialSetAtoms,
   METH_VARARGS, "Set the atoms prior to the first calculation."},
  {"get_potential_energy", (PyCFunction)PyAsap_PotentialGetPotentialEnergy,
   METH_VARARGS, "Calculate the potential energy."},
  {"get_potential_energies", (PyCFunction)PyAsap_PotentialGetPotentialEnergies,
   METH_VARARGS, "Calculate the potential energies of all atoms."},
  {"get_forces", (PyCFunction)PyAsap_PotentialGetForces,
   METH_VARARGS, "Calculate the potential energies of all atoms."},
  {"get_stress", (PyCFunction)PyAsap_PotentialGetStress,
   METH_VARARGS, "Calculate the potential energies of all atoms."},
  {"get_stresses", (PyCFunction)PyAsap_PotentialGetStresses,
   METH_VARARGS, "Calculate the potential energies of all atoms."},
  {"calculation_required", (PyCFunction)PyAsap_CalcReq,
   METH_VARARGS, "Check if a calculation is required."},
  {"get_neighborlist", (PyCFunction)PyAsap_PotentialGetNbList,
   METH_NOARGS,  "Return the neighbor list (if any)."},
  {"supports_parallel", (PyCFunction)PyAsap_Parallel,
   METH_NOARGS,  "Return True if the calculator supports parallel simulations."},
  {"get_cutoff", (PyCFunction)PyAsap_PotentialGetCutoff,
   METH_NOARGS, "Get the cutoff distance published by the calculator."},
  {"print_memory", (PyCFunction)PyAsap_PotPrintMemory,
   METH_NOARGS,  "Print an estimate of the memory usage."},
  {NULL}  // Sentinel
};

static PyMethodDef PyAsap_EMTMethods[] = {
  {"get_extra", (PyCFunction)PyAsap_PotentialGetExtra,
   METH_VARARGS, "Return extra information about the potential"},
  {"set_subtractE0", (PyCFunction)PyAsap_EMTSetSubtractE0,
   METH_VARARGS, "Set to false to move energy zero-point to infinitely separated atoms."},
  {"_use_imageatoms", (PyCFunction) PyAsap_PotentialUseImageAtoms,
   METH_NOARGS, PyAsap_PotentialUseImageAtoms_Docstring},
  {NULL}  // Sentinel
};

static PyMethodDef PyAsap_EMT2013Methods[] = {
  {"get_extra", (PyCFunction)PyAsap_PotentialGetExtra,
   METH_VARARGS, "Return extra information about the potential"},
  {NULL}  // Sentinel
};

namespace ASAPSPACE {
void InitPotentialType(PyTypeObject &type, bool mark)
{
  type.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE;
  type.tp_base = &PyAsap_PotentialType;
}
			 
int PyAsap_InitPotentialInterface(PyObject *module)
{
  // Initialize static string
  PyAsap_String_Energy = PyString_FromString("energy");
  PyAsap_String_Forces = PyString_FromString("forces");
  PyAsap_String_Stress = PyString_FromString("stress");

  // Init the potentials
  PyAsap_PotentialType.tp_new = PyType_GenericNew;
  PyAsap_PotentialType.tp_dealloc = PyAsap_PotentialDealloc;
  PyAsap_PotentialType.tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE;
  PyAsap_PotentialType.tp_methods = PyAsap_PotentialMethods;
  PyAsap_PotentialType.tp_repr = PyAsap_Representation<PyAsap_PotentialObject>;
  PyAsap_PotentialType.tp_init = (initproc) PyAsap_PotentialInit;
  PyAsap_PotentialType.tp_doc = Potential_Docstring;
  if (PyType_Ready(&PyAsap_PotentialType) < 0)
    return -1;

  InitPotentialType(PyAsap_EMTType);
  PyAsap_EMTType.tp_init = (initproc) PyAsap_EMTInit;
  PyAsap_EMTType.tp_doc = EMT_Docstring;
  PyAsap_EMTType.tp_methods = PyAsap_EMTMethods;
  if (PyType_Ready(&PyAsap_EMTType) < 0)
    return -1;
  Py_INCREF(&PyAsap_EMTType);
  PyModule_AddObject(module, "EMT", (PyObject *) &PyAsap_EMTType);

  InitPotentialType(PyAsap_EMT2013Type);
  PyAsap_EMT2013Type.tp_init = (initproc) PyAsap_EMT2013Init;
  PyAsap_EMT2013Type.tp_doc = EMT2013_Docstring;
  PyAsap_EMTType.tp_methods = PyAsap_EMT2013Methods;
  if (PyType_Ready(&PyAsap_EMT2013Type) < 0)
    return -1;
  Py_INCREF(&PyAsap_EMT2013Type);
  PyModule_AddObject(module, "EMT2013", (PyObject *) &PyAsap_EMT2013Type);

#ifdef MCEMT
  InitPotentialType(PyAsap_MonteCarloEMTType);
  PyAsap_MonteCarloEMTType.tp_init = (initproc) PyAsap_MonteCarloEMTInit;
  PyAsap_MonteCarloEMTType.tp_doc = MonteCarloEMT_Docstring;
  if (PyType_Ready(&PyAsap_MonteCarloEMTType) < 0)
    return -1;
  Py_INCREF(&PyAsap_MonteCarloEMTType);
  PyModule_AddObject(module, "MonteCarloEMT",
		     (PyObject *) &PyAsap_MonteCarloEMTType);
#endif
  
  InitPotentialType(PyAsap_RGLType);
  PyAsap_RGLType.tp_init = (initproc) PyAsap_RGLInit;
  PyAsap_RGLType.tp_doc = RGL_Docstring;
  if (PyType_Ready(&PyAsap_RGLType) < 0)
    return -1;
  Py_INCREF(&PyAsap_RGLType);
  PyModule_AddObject(module, "RGL", (PyObject *) &PyAsap_RGLType);

  InitPotentialType(PyAsap_LennardJonesType);
  PyAsap_LennardJonesType.tp_init = (initproc) PyAsap_LennardJonesInit;
  PyAsap_LennardJonesType.tp_doc = LennardJones_Docstring;
  if (PyType_Ready(&PyAsap_LennardJonesType) < 0)
    return -1;
  Py_INCREF(&PyAsap_LennardJonesType);
  PyModule_AddObject(module, "LennardJones",
		     (PyObject *) &PyAsap_LennardJonesType);

  InitPotentialType(PyAsap_MorseType);
  PyAsap_MorseType.tp_init = (initproc) PyAsap_MorseInit;
  PyAsap_MorseType.tp_doc = Morse_Docstring;
  if (PyType_Ready(&PyAsap_MorseType) < 0)
    return -1;
  Py_INCREF(&PyAsap_MorseType);
  PyModule_AddObject(module, "Morse",
		     (PyObject *) &PyAsap_MorseType);

  BrennerPotential::Initialize();
  InitPotentialType(PyAsap_BrennerType);
  PyAsap_BrennerType.tp_init = (initproc) PyAsap_BrennerInit;
  PyAsap_BrennerType.tp_doc = Brenner_Docstring;
  if (PyType_Ready(&PyAsap_BrennerType) < 0)
    return -1;
  Py_INCREF(&PyAsap_BrennerType);
  PyModule_AddObject(module, "BrennerPotential",
		     (PyObject *) &PyAsap_BrennerType);
return 0;
}

bool PyAsap_PotentialCheck(PyObject *obj)
{
  return PyObject_TypeCheck(obj, &PyAsap_PotentialType);
}

} // end namespace
