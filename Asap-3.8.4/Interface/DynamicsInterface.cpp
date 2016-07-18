// DynamicsInterface.cpp: Python interface to the dynamics objects.
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

#include "DynamicsInterface.h"
#include "VelocityVerlet.h"
#include "Langevin.h"
#include "ExceptionInterface.h"
#include "PythonConversions.h"
#include "Templates.h"
#include "Potential.h"

namespace ASAPSPACE {

static PyTypeObject PyAsap_VelocityVerletType = {
  PyObject_HEAD_INIT(NULL)
  0,
  "asapserial3.VelocityVerlet",
  sizeof(PyAsap_DynamicsObject),
  // The rest are initialized by name for reliability.
};
  
static PyTypeObject PyAsap_LangevinType = {
  PyObject_HEAD_INIT(NULL)
  0,
  "asapserial3.Langevin",
  sizeof(PyAsap_DynamicsObject),
  // The rest are initialized by name for reliability.
};

static char VelocityVerlet_Docstring[] = "ASAP-optimized Velocity Verlet dynamics object.\n";

static char Langevin_Docstring[] = "Asap-optimized Langevin dynamics object.\n";

// A few convenience macros

#define CHECK_DYN_UNINIT if (self->cobj != NULL) {		\
    PyErr_SetString(PyAsap_ErrorObject,				\
		    "Dynamics object already initialized.");	\
    return -1;							\
  }

#define CHECK_DYN_INIT if (self->cobj == NULL) {	  \
    PyErr_SetString(PyAsap_ErrorObject,			  \
		    "Dynamics object not initialized."); \
    return NULL;					  \
  }
  
static int PyAsap_VelocityVerletInit(PyAsap_DynamicsObject *self, 
    PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {"atoms", "calc", "timestep", NULL};
  
    PyObject *atoms;
    PyObject *calc;
    double timestep;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs,  "OOd:VelocityVerlet", kwlist,
            &atoms, &calc, &timestep))
        return -1;
    CHECK_DYN_UNINIT;
    Potential *asap_calc = NULL;
    if (PyAsap_PotentialCheck(calc))
    {
        // We can use the potential directly.
        asap_calc = ((PyAsap_PotentialObject *)calc)->cobj;
    }
    
    try
    {
        self->cobj = new VelocityVerlet(atoms, asap_calc, timestep);
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

static int PyAsap_LangevinInit(PyAsap_DynamicsObject *self,
    PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"atoms", "calc", "timestep", "sdpos", "sdmom", "c1", "c2", "fixcm", "seed", NULL};

  PyObject *atoms;
  PyObject *calc;
  double timestep;
  PyObject *sdpos;
  PyObject *sdmom;
  PyObject *c1;
  PyObject *c2;
  int fixcm;
  unsigned int seed;
  if (!PyArg_ParseTupleAndKeywords(args, kwargs,  "OOdO!O!O!O!iI:Langevin", kwlist,
          &atoms, &calc, &timestep, &PyString_Type, &sdpos, &PyString_Type, &sdmom,
	  &PyString_Type, &c1, &PyString_Type, &c2, &fixcm, &seed))
      return -1;
  CHECK_DYN_UNINIT;
  Potential *asap_calc = NULL;
  if (PyAsap_PotentialCheck(calc))
  {
      // We can use the potential directly.
      asap_calc = ((PyAsap_PotentialObject *)calc)->cobj;
  }

  try
  {
      self->cobj = new Langevin(atoms, asap_calc, timestep, sdpos, sdmom, c1,
				c2, (bool) fixcm, seed);
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

static PyObject *PyAsap_DynamicsRun(PyAsap_DynamicsObject *self,
    PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {"steps", "observers", "dyn", NULL};
    int steps;
    PyObject *observers;
    PyObject *dyn;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs,  "iOO:run", kwlist,
            &steps, &observers, &dyn))
        return NULL;
    CHECK_DYN_INIT;
    try
    {
        self->cobj->Run(steps, observers, dyn);
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

static PyObject *PyAsap_LgvSetScalarConstants(PyAsap_DynamicsObject *self,
    PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {"act0", "c3", "c4", "pmcor", "cnst", NULL};
    double act0;
    double c3;
    double c4;
    double pmcor;
    double cnst;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs,  "ddddd:set_scalar_constants", kwlist,
            &act0, &c3, &c4, &pmcor, &cnst))
        return NULL;
    CHECK_DYN_INIT;
    Langevin *cobj = dynamic_cast<Langevin *>(self->cobj);
    if (cobj == NULL)
      {
        PyErr_SetString(PyExc_TypeError, "Apparently not a Langevin object.");
        return NULL;
      }
    // Cannot throw (useful) exceptions.
    cobj->SetScalarConstants(act0, c3, c4, pmcor, cnst);
    Py_RETURN_NONE;
}

static PyObject *PyAsap_LgvSetVectorConstants(PyAsap_DynamicsObject *self,
    PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {"act0", "c3", "c4", "pmcor", "cnst", NULL};
    PyObject *act0;
    PyObject *c3;
    PyObject *c4;
    PyObject *pmcor;
    PyObject *cnst;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs,  "O!O!O!O!O!:set_vector_constants", kwlist,
        &PyString_Type, &act0, &PyString_Type, &c3, &PyString_Type, &c4,
        &PyString_Type, &pmcor, &PyString_Type, &cnst))
        return NULL;
    CHECK_DYN_INIT;
    Langevin *cobj = dynamic_cast<Langevin *>(self->cobj);
    if (cobj == NULL)
      {
        PyErr_SetString(PyExc_TypeError, "Apparently not a Langevin object.");
        return NULL;
      }
    try
    {
        cobj->SetVectorConstants(act0, c3, c4, pmcor, cnst);
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

static PyObject *PyAsap_LgvGetRandom(PyAsap_DynamicsObject *self,
                                     PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {"gaussian", NULL};
  int gaussian;
  if (!PyArg_ParseTupleAndKeywords(args, kwargs,  "i:get_random", kwlist,
      &gaussian))
      return NULL;
  CHECK_DYN_INIT;
  Langevin *cobj = dynamic_cast<Langevin *>(self->cobj);
  if (cobj == NULL)
    {
      PyErr_SetString(PyExc_TypeError, "Apparently not a Langevin object.");
      return NULL;
    }
  vector<Vec> x1;
  vector<Vec> x2;
  cobj->GetRandom(x1, x2, (bool) gaussian);
  return Py_BuildValue("NN",
      PyAsap_ArrayFromVectorVec(x1),
      PyAsap_ArrayFromVectorVec(x2));
}

static PyMethodDef PyAsap_VelocityVerletMethods[] = {
    {"run", (PyCFunction)PyAsap_DynamicsRun,
        METH_VARARGS|METH_KEYWORDS, "Run the dynamics."},
    {NULL}  // Sentinel
};

static PyMethodDef PyAsap_LangevinMethods[] = {
    {"run", (PyCFunction)PyAsap_DynamicsRun,
        METH_VARARGS|METH_KEYWORDS, "Run the dynamics."},
    {"set_scalar_constants", (PyCFunction)PyAsap_LgvSetScalarConstants,
        METH_VARARGS|METH_KEYWORDS, "Set constants (scalar version)."},
    {"set_vector_constants", (PyCFunction)PyAsap_LgvSetVectorConstants,
        METH_VARARGS|METH_KEYWORDS, "Set constants (scalar version)."},
    {"get_random", (PyCFunction)PyAsap_LgvGetRandom,
        METH_VARARGS|METH_KEYWORDS, "Get two sets of random numbers"},
    {NULL}  // Sentinel
};

static void InitDynamicsType(PyTypeObject &type)
{
  type.tp_new = PyType_GenericNew;
  type.tp_dealloc = PyAsap_Dealloc<PyAsap_DynamicsObject>;
  type.tp_flags = Py_TPFLAGS_DEFAULT;
  type.tp_repr = PyAsap_Representation<PyAsap_DynamicsObject>;
}
			 
int PyAsap_InitDynamicsInterface(PyObject *module)
{
    // Init the Velocity Verlet dynamics
    InitDynamicsType(PyAsap_VelocityVerletType);
    PyAsap_VelocityVerletType.tp_init = (initproc) PyAsap_VelocityVerletInit;
    PyAsap_VelocityVerletType.tp_doc = VelocityVerlet_Docstring;
    PyAsap_VelocityVerletType.tp_methods = PyAsap_VelocityVerletMethods;
    if (PyType_Ready(&PyAsap_VelocityVerletType) < 0)
        return -1;
    Py_INCREF(&PyAsap_VelocityVerletType);
    PyModule_AddObject(module, "VelocityVerlet", (PyObject *) 
        &PyAsap_VelocityVerletType);

    // Init the Langevin dynamics
    InitDynamicsType(PyAsap_LangevinType);
    PyAsap_LangevinType.tp_init = (initproc) PyAsap_LangevinInit;
    PyAsap_LangevinType.tp_doc = Langevin_Docstring;
    PyAsap_LangevinType.tp_methods = PyAsap_LangevinMethods;
    if (PyType_Ready(&PyAsap_LangevinType) < 0)
        return -1;
    Py_INCREF(&PyAsap_LangevinType);
    PyModule_AddObject(module, "Langevin", (PyObject *)
        &PyAsap_LangevinType);

    return 0;
}
        
} // end namespace
