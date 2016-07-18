// -*- C++ -*-
//
// Copyright (C) 2008-2011 Jakob Schiotz and Center for Individual
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


#include "AsapPython.h"
#include <structmember.h>
#include <mpi.h>
#include "mpimodule.h"

#if ASAPDEBUG
#include <stdio.h>
#define DEBUGPRINT fprintf(stderr, "%s:%d  %s\n", __FILE__, __LINE__, __FUNCTION__)
#else
#define DEBUGPRINT
#endif

// #define LONGP(a) ((long*)((a)->data))
#define GPAW_MALLOC(T, n) ((T*)malloc((n) * sizeof(T)))

// Check that array is well-behaved and contains data that can be sent.
#define CHK_ARRAY(a) if ((a) == NULL                    		\
			 || !PyArray_ISCARRAY(a) || !PyArray_ISNUMBER(a)) { \
    PyErr_SetString(PyExc_TypeError,					\
		    "Not a proper NumPy array for MPI communication."); \
    return NULL; }

#define CHK_ARRAY_2(a) if ((a) == NULL || !PyArray_Check(a)             \
                         || !PyArray_ISCARRAY(a) || !PyArray_ISNUMBER(a)) { \
    PyErr_SetString(PyExc_TypeError,                                    \
                    "Not a proper NumPy array for MPI communication."); \
    return NULL; }

// Check that two arrays have the same type, and the size of the
// second is a given multiple of the size of the first
#define CHK_ARRAYS(a,b,n)						\
  if ((PyArray_TYPE(a) != PyArray_TYPE(b))				\
      || (PyArray_SIZE(b) != PyArray_SIZE(a) * n)) {			\
    PyErr_SetString(PyExc_ValueError,					\
		    "Incompatible array types or sizes.");		\
      return NULL; }

// Check that a processor number is valid
#define CHK_PROC(n) if (n < 0 || n >= self->size) {\
    PyErr_SetString(PyExc_ValueError, "Invalid processor number.");	\
    return NULL; }

// Check that a processor number is valid or is -1
#define CHK_PROC_DEF(n) if (n < -1 || n >= self->size) {\
    PyErr_SetString(PyExc_ValueError, "Invalid processor number.");	\
    return NULL; }

// Check that a processor number is valid and is not this processor
#define CHK_OTHER_PROC(n) if (n < 0 || n >= self->size || n == self->rank) { \
    PyErr_SetString(PyExc_ValueError, "Invalid processor number.");	\
    return NULL; }

// MPI request object, so we can store a reference to the buffer,
// preventing its early deallocation.
typedef struct {
  PyObject_HEAD
  MPI_Request rq;
  PyArrayObject *buffer;
  int status;
} mpi_request;

static PyObject *mpi_request_wait(mpi_request *self, PyObject *noargs)
{
  DEBUGPRINT;
  if (self->status == 0)
    {
      // Calling wait multiple times is allowed but meaningless (as in the MPI standard)
      Py_RETURN_NONE;
    }
  MPI_Wait(&(self->rq), MPI_STATUS_IGNORE);
  CHECKREF(self->buffer);
  Py_DECREF(self->buffer);
  self->status = 0;
  DEBUGPRINT;
  Py_RETURN_NONE;
}

static PyObject *mpi_request_test(mpi_request *self, PyObject *noargs)
{
  DEBUGPRINT;
  if (self->status == 0)
    {
      Py_RETURN_TRUE;  // Already completed
    }
  int flag;
  MPI_Test(&(self->rq), &flag, MPI_STATUS_IGNORE);
  if (flag)
    {
      CHECKREF(self->buffer);
      Py_DECREF(self->buffer);
      self->status = 0;
      DEBUGPRINT;
      Py_RETURN_TRUE;
    }
  else
    {
      Py_RETURN_FALSE;
    }
}

static void mpi_request_dealloc(mpi_request *self)
{
  DEBUGPRINT;
  if (self->status)
    {
      DEBUGPRINT;
      mpi_request_wait(self, NULL);
    }
  PyObject_Del(self);
  DEBUGPRINT;
}

static PyMemberDef mpi_request_members[] = {
    {"status", T_INT, offsetof(mpi_request, status), READONLY,
        "status of the request, non-zero if communication is pending."},
    {NULL}
};

static PyMethodDef mpi_request_methods[] = {
    {"wait", (PyCFunction) mpi_request_wait, METH_NOARGS,
        "Wait for the communication to complete."
    },
    {"test", (PyCFunction) mpi_request_test, METH_NOARGS,
        "Test if the communication has completed."
    },
    {NULL}
};

PyTypeObject mpi_request_type = {
    PyObject_HEAD_INIT(&PyType_Type)
    0,                         /*ob_size*/
    "MPI_Request",             /*tp_name*/
    sizeof(mpi_request),             /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)mpi_request_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "MPI request object",           /* tp_doc */
    0,                   /* tp_traverse */
    0,                   /* tp_clear */
    0,                   /* tp_richcompare */
    0,                   /* tp_weaklistoffset */
    0,                   /* tp_iter */
    0,                   /* tp_iternext */
    mpi_request_methods,             /* tp_methods */
    mpi_request_members,
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    0,      /* tp_init */
    0,                         /* tp_alloc */
    0,                 /* tp_new */

};


static mpi_request *PyAsap_NewMPIRequest()
{
  mpi_request *self;
  DEBUGPRINT;

  self = PyObject_NEW(mpi_request, &mpi_request_type);
  if (self == NULL) return NULL;
  self->buffer = NULL;
  self->status = 1;  // Active
  DEBUGPRINT;
  return self;
}

static void mpi_dealloc(MPIObject *obj)
{
  if (obj->comm != MPI_COMM_WORLD)
    {
      MPI_Comm_free(&(obj->comm));
      CHECKREF(obj->parent);
      Py_XDECREF(obj->parent);
    }
  PyObject_DEL(obj);
}

static PyObject * mpi_receive(MPIObject *self, PyObject *args, PyObject *kwargs)
{
  PyArrayObject* a;
  int src = -1;
  int tag = 123;
  int block = 1;
  static char *kwlist[] = {"a", "src", "tag", "block", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!i|ii:receive", kwlist,
				   &PyArray_Type, &a, &src, &tag, &block))
    return NULL;
  CHK_ARRAY(a);
  if (src != -1)
    {
      CHK_OTHER_PROC(src);
    }
  else
    src = MPI_ANY_SOURCE;
  int n = PyArray_DESCR(a)->elsize;
  int d;
  for (d = 0; d < PyArray_NDIM(a); d++)
    n *= PyArray_DIM(a, d);
  if (block)
    {
      MPI_Status status;
      MPI_Recv(PyArray_BYTES(a), n, MPI_BYTE, src, tag, self->comm,
	       &status);
      int real_source = status.MPI_SOURCE;
      return PyInt_FromLong((long) real_source);
    }
  else
    {
      mpi_request *req = PyAsap_NewMPIRequest();
      if (req == NULL) return NULL;
      req->buffer = a;
      Py_INCREF(req->buffer);
      MPI_Irecv(PyArray_BYTES(a), n, MPI_BYTE, src, tag, self->comm, &(req->rq));
      return (PyObject *) req;
    }
}

static PyObject * mpi_probe(MPIObject *self, PyObject *args, PyObject *kwargs)
{
  int src = -1;
  int tag = -1;
  int block = 0;
  static char *kwlist[] = {"src", "tag", "block", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|iii:probe", kwlist,
				   &src, &tag, &block))
    return NULL;
  
  if (src != -1)
    {
      CHK_OTHER_PROC(src);
    }
  else
    src = MPI_ANY_SOURCE;
  if (tag == -1)
    tag = MPI_ANY_TAG;

  MPI_Status status;
  int flag;
  if (block)
    {
      MPI_Probe(src, tag, self->comm, &status);
      flag = 1;
    }
  else
    {
      MPI_Iprobe(src, tag, self->comm, &flag, &status);
    }
  if (flag)
    {
      src = status.MPI_SOURCE;
      tag = status.MPI_TAG;
      int count;
      MPI_Get_count(&status, MPI_BYTE, &count);
      return Py_BuildValue("iii", src, tag, count);
    }
  else
    Py_RETURN_NONE;
}

static PyObject * mpi_send(MPIObject *self, PyObject *args, PyObject *kwargs)
{
  PyArrayObject* a;
  int dest;
  int tag = 123;
  int block = 1;
  static char *kwlist[] = {"a", "dest", "tag", "block", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!i|ii:send", kwlist,
				   &PyArray_Type, &a, &dest, &tag, &block))
    return NULL;
  CHK_ARRAY(a);
  CHK_OTHER_PROC(dest);
  int n = PyArray_DESCR(a)->elsize;
  int d;
  for (d = 0; d < PyArray_NDIM(a); d++)
    n *= PyArray_DIM(a,d);
  if (block)
    {
      MPI_Send(PyArray_BYTES(a), n, MPI_BYTE, dest, tag, self->comm);
      Py_RETURN_NONE;
    }
  else
    {
      mpi_request *req = PyAsap_NewMPIRequest();
      if (req == NULL) return NULL;
      req->buffer = a;
      Py_INCREF(a);
      MPI_Isend(PyArray_BYTES(a), n, MPI_BYTE, dest, tag, self->comm,
		&(req->rq));
      return (PyObject *) req;
    }
}


static PyObject * mpi_sendreceive(MPIObject *self, PyObject *args,
                                  PyObject *kwargs)
{
  PyArrayObject* a;
  PyArrayObject* b;
  int dest, src;
  int sendtag = 123;
  int recvtag = 123;
  static char *kwlist[] = {"a", "dest", "b", "src", "sendtag", "recvtag",
      NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!iO!i|ii:sendreceive",
      kwlist,
      &PyArray_Type, &a, &dest, &PyArray_Type, &b, &src, &sendtag, &recvtag))
    return NULL;
  CHK_ARRAY(a);
  CHK_OTHER_PROC(dest);
  CHK_ARRAY(b);
  CHK_OTHER_PROC(src);
  int nsend = PyArray_DESCR(a)->elsize;
  int d;
  for (d = 0; d < PyArray_NDIM(a); d++)
    nsend *= PyArray_DIM(a,d);
  int nrecv = PyArray_DESCR(b)->elsize;
  for (d = 0; d < PyArray_NDIM(b); d++)
    nrecv *= PyArray_DIM(b,d);
  MPI_Sendrecv(PyArray_BYTES(a), nsend, MPI_BYTE, dest, sendtag,
               PyArray_BYTES(b), nrecv, MPI_BYTE, src, recvtag,
               self->comm, MPI_STATUS_IGNORE);
  Py_RETURN_NONE;
}

static PyObject * mpi_name(MPIObject *self, PyObject *noargs)
{
  char name[MPI_MAX_PROCESSOR_NAME];
  int resultlen;
  MPI_Get_processor_name(name, &resultlen);
  return Py_BuildValue("s#", name, resultlen);
}

static PyObject * mpi_abort(MPIObject *self, PyObject *args)
{
  int errcode;
  if (!PyArg_ParseTuple(args, "i:abort", &errcode))
    return NULL;
  MPI_Abort(self->comm, errcode);
  Py_RETURN_NONE;
}

static PyObject * mpi_barrier(MPIObject *self)
{
  MPI_Barrier(self->comm);
  Py_RETURN_NONE;
}

static PyObject * mpi_wait(MPIObject *self, PyObject *args)
{
  mpi_request* s;
  if (!PyArg_ParseTuple(args, "O!:wait", &mpi_request_type, &s))
    return NULL;
  return mpi_request_wait(s, NULL);
}

static PyObject * mpi_test(MPIObject *self, PyObject *args)
{
  mpi_request* s;
  if (!PyArg_ParseTuple(args, "O!:wait", &mpi_request_type, &s))
        return NULL;
  return mpi_request_test(s, NULL);
}

static PyObject * mpi_waitall(MPIObject *self, PyObject *requests)
{
  int n;   // Number of requests
  MPI_Request *rqs = NULL;
  if (!PySequence_Check(requests))
    {
      PyErr_SetString(PyExc_TypeError, "mpi.waitall: argument must be a sequence");
      return NULL;
    }
  // Extract the request objects
  n = PySequence_Size(requests);
  assert(n >= 0);  // This cannot fail.
  rqs = GPAW_MALLOC(MPI_Request, n);
  int i;
  for (i = 0; i < n; i++)
    {
      PyObject *o = PySequence_GetItem(requests, i);
      if (o == NULL)
        return NULL;
      if (o->ob_type != &mpi_request_type)
        {
          Py_DECREF(o);
          free(rqs);
          PyErr_SetString(PyExc_TypeError, "mpi.waitall: argument must be a sequence of MPI requests");
          return NULL;
        }
      mpi_request *s = (mpi_request *)o;
      rqs[i] = s->rq;
      Py_DECREF(o);
    }
  // Do the actual wait.
  MPI_Waitall(n, rqs, MPI_STATUSES_IGNORE);

  // Release the buffers used by the MPI communication
  for (i = 0; i < n; i++)
   {
     mpi_request *o = (mpi_request *) PySequence_GetItem(requests, i);
     if (o->status)
     {
       assert(o->buffer != NULL);
       Py_DECREF(o->buffer);
     }
     o->status = 0;
     Py_DECREF(o);
   }
  // Release internal data and return.
  free(rqs);
  Py_RETURN_NONE;
}

static PyObject * mpi_testall(MPIObject *self, PyObject *requests)
{
  int n;   // Number of requests
  MPI_Request *rqs = NULL;
  int flag = 0;
  DEBUGPRINT;
  if (!PySequence_Check(requests))
    {
      PyErr_SetString(PyExc_TypeError, "mpi.testall: argument must be a sequence");
      return NULL;
    }
  // Extract the request objects
  n = PySequence_Size(requests);
  assert(n >= 0);  // This cannot fail.
  rqs = GPAW_MALLOC(MPI_Request, n);
  assert(rqs != NULL);
  DEBUGPRINT;
  int i;
  for (i = 0; i < n; i++)
    {
      PyObject *o = PySequence_GetItem(requests, i);
      if (o == NULL)
        return NULL;
      if (o->ob_type != &mpi_request_type)
        {
          Py_DECREF(o);
          free(rqs);
          PyErr_SetString(PyExc_TypeError, "mpi.testall: argument must be a sequence of MPI requests");
          return NULL;
        }
      mpi_request *s = (mpi_request *)o;
      rqs[i] = s->rq;
      Py_DECREF(o);
    }
  // Do the actual test.
  DEBUGPRINT;
  MPI_Testall(n, rqs, &flag, MPI_STATUSES_IGNORE);
  // Unlike MPI_Test, if flag outcome is non-zero, MPI_Testall will deallocate
  // all requests which were allocated by nonblocking communication calls, so
  // we must free these buffers. Otherwise, none of the requests are modified.
  DEBUGPRINT;
  if (flag != 0)
    {
      // Release the buffers used by the MPI communication
    DEBUGPRINT;
      for (i = 0; i < n; i++)
      {
        mpi_request *o = (mpi_request *) PySequence_GetItem(requests, i);
        if (o->status)
        {
          assert(o->buffer != NULL);
          Py_DECREF(o->buffer);
        }
        o->status = 0;
        Py_DECREF(o);
      }
    }
  // Release internal data and return.
  DEBUGPRINT;
  free(rqs);
  DEBUGPRINT;
  return Py_BuildValue("i", flag);
}


static MPI_Datatype get_mpi_datatype(PyArrayObject *a)
{
  int n = PyArray_DESCR(a)->elsize;
  if (PyArray_ISCOMPLEX(a))
    n = n/2;

  switch(PyArray_TYPE(a))
    {
      // Floating point numbers including complex numbers
    case NPY_DOUBLE:
    case NPY_CDOUBLE:
      assert(sizeof(double) == n);
      return MPI_DOUBLE;
    case NPY_FLOAT:
    case NPY_CFLOAT:
      assert(sizeof(float) == n);
      return MPI_FLOAT;
    case NPY_LONGDOUBLE:
    case NPY_CLONGDOUBLE:
       assert(sizeof(long double) == n);
      return MPI_LONG_DOUBLE;
      // Signed integer types
    case NPY_BYTE:
      assert(sizeof(char) == n);
      return MPI_CHAR;
    case NPY_SHORT:
      assert(sizeof(short) == n);
      return MPI_SHORT;
    case NPY_INT:
      assert(sizeof(int) == n);
      return MPI_INT;
    case NPY_LONG:
      assert(sizeof(long) == n);
      return MPI_LONG;
      // Unsigned integer types
    case NPY_BOOL:
    case NPY_UBYTE:
      assert(sizeof(unsigned char) == n);
      return MPI_UNSIGNED_CHAR;
    case NPY_USHORT:
      assert(sizeof(unsigned short) == n);
      return MPI_UNSIGNED_SHORT;
    case NPY_UINT:
      assert(sizeof(unsigned) == n);
      return MPI_UNSIGNED;
    case NPY_ULONG:
      assert(sizeof(unsigned long) == n);
      return MPI_UNSIGNED_LONG;
    }
  // If we reach this point none of the cases worked out.
  PyErr_SetString(PyExc_ValueError, "Cannot communicate data of this type.");
  return 0;
}

static PyObject * mpi_reduce(MPIObject *self, PyObject *args, PyObject *kwargs,
			     MPI_Op operation, int allowcomplex)
{
  PyObject* obj;
  int root = -1;
  static char *kwlist[] = {"a", "root", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|i:reduce", kwlist,
				   &obj, &root))
    return NULL;
  CHK_PROC_DEF(root);
  if (PyFloat_Check(obj))
    {
      double din = PyFloat_AS_DOUBLE(obj);
      double dout;
      if (root == -1)
        MPI_Allreduce(&din, &dout, 1, MPI_DOUBLE, operation, self->comm);
      else
        MPI_Reduce(&din, &dout, 1, MPI_DOUBLE, operation, root, self->comm);
      return PyFloat_FromDouble(dout);
    }
  if (PyInt_Check(obj))
    {
      long din = PyInt_AS_LONG(obj);
      long dout;
      if (root == -1)
        MPI_Allreduce(&din, &dout, 1, MPI_LONG, operation, self->comm);
      else
        MPI_Reduce(&din, &dout, 1, MPI_LONG, operation, root, self->comm);
      return PyInt_FromLong(dout);
    }
  else if (PyComplex_Check(obj) && allowcomplex)
    {
      double din[2];
      double dout[2];
      din[0] = PyComplex_RealAsDouble(obj);
      din[1] = PyComplex_ImagAsDouble(obj);
      if (root == -1)
        MPI_Allreduce(&din, &dout, 2, MPI_DOUBLE, MPI_SUM, self->comm);
      else
        MPI_Reduce(&din, &dout, 2, MPI_DOUBLE, MPI_SUM, root, self->comm);
      return PyComplex_FromDoubles(dout[0], dout[1]);
    }
  else if (PyComplex_Check(obj))
    {
      PyErr_SetString(PyExc_ValueError,
		      "Operation not allowed on complex numbers");
      return NULL;
    }   
  else   // It should be an array
    {
      PyArrayObject *arr = (PyArrayObject *) obj;
      int n;
      int elemsize;
      MPI_Datatype datatype;
      CHK_ARRAY_2(arr);
      datatype = get_mpi_datatype(arr);
      if (datatype == 0)
	return NULL;
      n = PyArray_SIZE(arr);
      elemsize = PyArray_DESCR(arr)->elsize;
      if (PyArray_ISCOMPLEX(arr))
	{
	  if (allowcomplex)
	    {
	      n *= 2;
	      elemsize /= 2;
	    }
	  else
	    {
	      PyErr_SetString(PyExc_ValueError,
			      "Operation not allowed on complex numbers");
	      return NULL;
	    }
	}
      if (root == -1)
        {
          char* b = GPAW_MALLOC(char, n * elemsize);
          // XXX Use MPI_IN_PLACE!!
          MPI_Allreduce(PyArray_BYTES(arr), b, n, datatype, operation,
			self->comm);
	  assert(PyArray_NBYTES(arr) == n * elemsize);
          memcpy(PyArray_BYTES(arr), b, n * elemsize);
          free(b);
        }
      else
        {
          char* b = 0;
          int rank;
          MPI_Comm_rank(self->comm, &rank);
#ifdef GPAW_BGP
          b = GPAW_MALLOC(char, n * elemsize); // bug on BGP
#else
          if (rank == root)
               b = GPAW_MALLOC(char, n * elemsize);
#endif
          // XXX Use MPI_IN_PLACE!!
          MPI_Reduce(PyArray_BYTES(arr), b, n, datatype, operation, root,
		     self->comm);
          if (rank == root)
            {
	      assert(PyArray_NBYTES(arr) == n * elemsize);
              memcpy(PyArray_BYTES(arr), b, n * elemsize);
            }
#ifdef GPAW_BGP
          free(b); // bug on BGP
#else
          if (rank == root)
               free(b);
#endif
        }
      Py_RETURN_NONE;
    }
}

static PyObject * mpi_sum(MPIObject *self, PyObject *args, PyObject *kwargs)
{
  return mpi_reduce(self, args, kwargs, MPI_SUM, 1);
}

static PyObject * mpi_product(MPIObject *self, PyObject *args, PyObject *kwargs)
{
  // No complex numbers as that would give separate products of
  // real and imaginary parts.
  return mpi_reduce(self, args, kwargs,  MPI_PROD, 0); 
}

static PyObject * mpi_max(MPIObject *self, PyObject *args, PyObject *kwargs)
{
  return mpi_reduce(self, args,  kwargs, MPI_MAX, 0);
}

static PyObject * mpi_min(MPIObject *self, PyObject *args, PyObject *kwargs)
{
  return mpi_reduce(self, args,  kwargs, MPI_MIN, 0);
}

static PyObject * mpi_scatter(MPIObject *self, PyObject *args)
{
  PyArrayObject* sendobj;
  PyArrayObject* recvobj;
  int root;
  if (!PyArg_ParseTuple(args, "O!O!i:scatter", &PyArray_Type, &sendobj,
                        &PyArray_Type, &recvobj, &root))
    return NULL;
  CHK_ARRAY(sendobj);
  CHK_ARRAY(recvobj);
  CHK_PROC(root);
  CHK_ARRAYS(recvobj, sendobj, self->size); // size(send) = size(recv)*Ncpu
  int n = PyArray_DESCR(recvobj)->elsize;
  int d;
  for (d = 0; d < PyArray_NDIM(recvobj); d++)
    n *= PyArray_DIM(recvobj,d);
  MPI_Scatter(PyArray_BYTES(sendobj), n, MPI_BYTE, PyArray_BYTES(recvobj),
	      n, MPI_BYTE, root, self->comm);
  Py_RETURN_NONE;
}



static PyObject * mpi_allgather(MPIObject *self, PyObject *args)
{
  PyArrayObject* a;
  PyArrayObject* b;
  if (!PyArg_ParseTuple(args, "O!O!:allgather",
                        &PyArray_Type, &a, &PyArray_Type, &b))
    return NULL;
  CHK_ARRAY(a);
  CHK_ARRAY(b);
  CHK_ARRAYS(a, b, self->size);
  int n = PyArray_DESCR(a)->elsize;
  int d;
  for (d = 0; d < PyArray_NDIM(a); d++)
    n *= PyArray_DIM(a,d);
  // What about endianness???? 
  MPI_Allgather(PyArray_BYTES(a), n, MPI_BYTE, PyArray_BYTES(b), n,
		MPI_BYTE, self->comm);
  Py_RETURN_NONE;
}

static PyObject * mpi_gather(MPIObject *self, PyObject *args)
{
  PyArrayObject* a;
  int root;
  PyArrayObject* b = 0;
  if (!PyArg_ParseTuple(args, "O!i|O!", &PyArray_Type, &a, &root,
                        &PyArray_Type, &b))
    return NULL;
  CHK_ARRAY(a);
  CHK_PROC(root);
  if (root == self->rank)
    {
      CHK_ARRAY(b);
      CHK_ARRAYS(a, b, self->size);
    }  // Else check for b=None or NULL.  But we could also just ignore b.
  int n = PyArray_DESCR(a)->elsize;
  int d;
  for (d = 0; d < PyArray_NDIM(a); d++)
    n *= PyArray_DIM(a,d);
  if (root != self->rank)  // What about endianness????
    MPI_Gather(PyArray_BYTES(a), n, MPI_BYTE, 0, n, MPI_BYTE, root, self->comm);
  else
    MPI_Gather(PyArray_BYTES(a), n, MPI_BYTE, PyArray_BYTES(b), n, MPI_BYTE, root, self->comm);
  Py_RETURN_NONE;
}

static PyObject * mpi_broadcast(MPIObject *self, PyObject *args)
{
  PyArrayObject* buf;
  int root;
  if (!PyArg_ParseTuple(args, "O!i:broadcast", &PyArray_Type, &buf, &root))
    return NULL;
  CHK_ARRAY(buf);
  CHK_PROC(root);
  int n = PyArray_DESCR(buf)->elsize;
  int d;
  for (d = 0; d < PyArray_NDIM(buf); d++)
    n *= PyArray_DIM(buf,d);
  MPI_Bcast(PyArray_BYTES(buf), n, MPI_BYTE, root, self->comm);
  Py_RETURN_NONE;
}

// Forward declaration of MPI_Communicator because it needs MPIType
// that needs MPI_getattr that needs MPI_Methods that need
// MPI_Communicator that need ...
static PyObject * MPICommunicator(MPIObject *self, PyObject *args);

static PyMethodDef mpi_methods[] = {
    {"receive",          (PyCFunction)mpi_receive,
     METH_VARARGS|METH_KEYWORDS,
     "receive(a, src, tag=123, block=1) receives array a from src."},
    {"send",             (PyCFunction)mpi_send,
     METH_VARARGS|METH_KEYWORDS,
     "send(a, dest, tag=123, block=1) sends array a to dest."},
    {"sendreceive",          (PyCFunction)mpi_sendreceive,
     METH_VARARGS|METH_KEYWORDS,
     "sendreceive(a, dest, b, src, sendtag=123, recvtag=123) sends an array a to dest and receives an array b from src."},
    {"probe",            (PyCFunction)mpi_probe, 
     METH_VARARGS|METH_KEYWORDS, 0},
    {"abort",            (PyCFunction)mpi_abort,        METH_VARARGS,
     "abort(errcode) aborts all MPI tasks."},
    {"name",             (PyCFunction)mpi_name,         METH_NOARGS,
     "name() returns the name of the processor node."},
    {"barrier",          (PyCFunction)mpi_barrier,      METH_VARARGS,
     "barrier() synchronizes all MPI tasks"},
    {"test",             (PyCFunction)mpi_test,         METH_VARARGS,
     "test(request) tests if a nonblocking communication is complete."},
    {"testall",          (PyCFunction)mpi_testall,      METH_O,
     "testall(list_of_rqs) tests if multiple nonblocking communications are complete."},
    {"wait",             (PyCFunction)mpi_wait,         METH_VARARGS,
     "wait(request) waits for a nonblocking communication to complete."},
    {"waitall",          (PyCFunction)mpi_waitall,      METH_O,
     "waitall(list_of_rqs) waits for multiple nonblocking communications to complete."},
    {"sum",              (PyCFunction)mpi_sum,
     METH_VARARGS|METH_KEYWORDS,
     "sum(a, root=-1) sums arrays or scalars, result on all tasks unless root is given."},
    {"product",          (PyCFunction)mpi_product,
     METH_VARARGS|METH_KEYWORDS,
     "product(a, root=-1) multiplies arrays or scalars, result on all tasks unless root is given."},
    {"max",              (PyCFunction)mpi_max,
     METH_VARARGS|METH_KEYWORDS,
     "max(a, root=-1) maximum of arrays or scalars, result on all tasks unless root is given."},
    {"min",              (PyCFunction)mpi_min,
     METH_VARARGS|METH_KEYWORDS,
     "min(a, root=-1) minimum of arrays or scalars, result on all tasks unless root is given."},
    {"scatter",          (PyCFunction)mpi_scatter,      METH_VARARGS,
     "scatter(src, target, root) distributes array from root task."},
    {"gather",           (PyCFunction)mpi_gather,       METH_VARARGS,
     "gather(src, root, target=None) gathers data from all tasks on root task."},
    {"all_gather",       (PyCFunction)mpi_allgather,    METH_VARARGS,
     "all_gather(src, target) gathers data from all tasks on all tasks."},
    {"broadcast",        (PyCFunction)mpi_broadcast,    METH_VARARGS, 0},
    {"new_communicator", (PyCFunction)MPICommunicator,  METH_VARARGS,
     "broadcast(buffer, root) Broadcast data in-place from root task."},
    {0, 0, 0, 0}
};

static PyMemberDef mpi_members[] = {
  {"size", T_INT, offsetof(MPIObject, size), 0, "Number of processors"},
  {"rank", T_INT, offsetof(MPIObject, rank), 0, "Number of this processor"},
  {0, 0, 0, 0, 0}  /* Sentinel */
};

// __new__
static PyObject *NewMPIObject(PyTypeObject* type, PyObject *args, PyObject *kwds)
{
  static char *kwlist[] = {NULL};
  MPIObject *self;
  
  if (! PyArg_ParseTupleAndKeywords(args, kwds, "", kwlist))
    return NULL;

  self = (MPIObject *) type->tp_alloc(type, 0);
  if (self == NULL)
    return NULL;
  MPI_Comm_size(MPI_COMM_WORLD, &(self->size));
  MPI_Comm_rank(MPI_COMM_WORLD, &(self->rank));
  self->comm = MPI_COMM_WORLD;
  self->parent = NULL;
  
  return (PyObject *) self;
}

// __init__
static int InitMPIObject(MPIObject* self, PyObject *args, PyObject *kwds)
{
  static char *kwlist[] = {NULL};

  if (! PyArg_ParseTupleAndKeywords(args, kwds, "", kwlist))
    return -1;

  return 0;
}

PyTypeObject MPIType = {
  PyObject_HEAD_INIT(&PyType_Type)
  0,                         /*ob_size*/
  "MPI",             /*tp_name*/
  sizeof(MPIObject),             /*tp_basicsize*/
  0,                         /*tp_itemsize*/
  (destructor)mpi_dealloc, /*tp_dealloc*/
  0,                         /*tp_print*/
  //  mpi_get_attr,                         /*tp_getattr*/
  0,                         /*tp_getattr*/
  0,                         /*tp_setattr*/
  0,                         /*tp_compare*/
  0,                         /*tp_repr*/
  0,                         /*tp_as_number*/
  0,                         /*tp_as_sequence*/
  0,                         /*tp_as_mapping*/
  0,                         /*tp_hash */
  0,                         /*tp_call*/
  0,                         /*tp_str*/
  0,                         /*tp_getattro*/
  0,                         /*tp_setattro*/
  0,                         /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
  "MPI object",           /* tp_doc */
  0,                   /* tp_traverse */
  0,                   /* tp_clear */
  0,                   /* tp_richcompare */
  0,                   /* tp_weaklistoffset */
  0,                   /* tp_iter */
  0,                   /* tp_iternext */
  mpi_methods,             /* tp_methods */
  mpi_members,
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)InitMPIObject,      /* tp_init */
    0,                         /* tp_alloc */
    NewMPIObject,                 /* tp_new */

};

static PyObject * MPICommunicator(MPIObject *self, PyObject *args)
{
  PyArrayObject* ranks;
  if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &ranks))
    return NULL;
  MPI_Group group;
  MPI_Comm_group(self->comm, &group);
  int n = PyArray_DIMS(ranks)[0];
  MPI_Group newgroup;
  // Stupid hack; MPI_Group_incl wants a int argument;
  // numpy arrays are long (might be different from ints)
  // More clever ways are welcomed...
  int* ranks_int = GPAW_MALLOC(int, n);
  long* ranks_long = PyArray_DATA(ranks);
  int i;
  for (i = 0; i < n ; i++ )
    ranks_int[i]=ranks_long[i];
  MPI_Group_incl(group, n, ranks_int, &newgroup);
  free(ranks_int);
  MPI_Comm comm;
  MPI_Comm_create(self->comm, newgroup, &comm);
  MPI_Group_free(&newgroup);
  MPI_Group_free(&group);
  if (comm == MPI_COMM_NULL)
    {
      Py_RETURN_NONE;
    }
  else
    {
      MPIObject *obj = PyObject_NEW(MPIObject, &MPIType);
      if (obj == NULL)
        return NULL;
      MPI_Comm_size(comm, &(obj->size));
      MPI_Comm_rank(comm, &(obj->rank));
      obj->comm = comm;
      // Make sure that MPI_COMM_WORLD is kept alive til the end (we
      // don't want MPI_Finalize to be called before MPI_Comm_free):
      Py_INCREF(self);
      obj->parent = (PyObject*)self;
      return (PyObject*)obj;
    }
}

int PyAsap_InitMpiInterface(PyObject *module)
{
  if (PyType_Ready(&MPIType) < 0)
    return -1;
  if (PyType_Ready(&mpi_request_type) < 0)
    return -1;

  Py_INCREF(&MPIType);
  Py_INCREF(&mpi_request_type);
  PyModule_AddObject(module, "Communicator", (PyObject *)&MPIType);

  return 0;
}
  
