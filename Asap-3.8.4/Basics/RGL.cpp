//#define ASAPDEBUG

#include "RGL.h"
#include "Atoms.h"
#include "NeighborCellLocator.h"
#include "NeighborList.h"
#include "Exception.h"
#include "Timing.h"
#include "mass.h"
#include "Debug.h"
#include <math.h>
#include <vector>
#include <iostream>
#include <algorithm>
using std::vector;
using std::cerr;
using std::endl;
using std::flush;

#undef INTERLEAVETHREADS

#if 1
#define VERB(x) if (verbose == 1) cerr << x
#else
#define VERB(x)
#endif

const static int stresscomp[3][3] = {{0, 5, 4}, {5, 1, 3}, {4, 3, 2}};

RGL::RGL (PyObject *self,
          const vector<int> &elements,
          const TinyMatrix<double> &p,
          const TinyMatrix<double> &q,
          const TinyMatrix<double> &A,
          const TinyMatrix<double> &qsi2,
          const TinyMatrix<double> &r0,
          const TinyMatrix<double> &p3,
          const TinyMatrix<double> &p4,
          const TinyMatrix<double> &p5,
          const TinyMatrix<double> &q3,
          const TinyMatrix<double> &q4,
          const TinyMatrix<double> &q5,
          const double &rcs,
          const double &rce) : Potential(self)
{
  DEBUGPRINT;
  CONSTRUCTOR;

  atoms = NULL;
  ghostatoms = false;
  neighborlist = NULL;
  neighborlist_obj = NULL;

  nAtoms = 0;
  nSize = 0;
  nElements = elements.size();
  driftfactor = 0.05;
  counters.energies = counters.forces = 0;

  param.p = p;
  param.q = q;
  param.A = A;
  param.qsi2 = qsi2;
  param.r0 = r0;
  param.p3 = p3;
  param.p4 = p4;
  param.p5 = p5;
  param.q3 = q3;
  param.q4 = q4;
  param.q5 = q5;

  this->cutoff.start = rcs;
  this->cutoff.end = rce;

  zmap.resize(92);
  for (int i = 0; i < 92; i++)
    zmap[i] = -1;
  for (int i = 0; i < nElements; i++)
    zmap[elements[i]] = i;
}

RGL::~RGL() {
  DESTRUCTOR;
  Py_XDECREF(neighborlist_obj);
  if (atoms != NULL)
    AsapAtoms_DECREF(atoms);
}

void RGL::SetAtoms(PyObject *pyatoms, Atoms* accessobj /* = NULL */) {
  DEBUGPRINT;

  // SetAtoms should only do anything the first time it is called.
  if (atoms != NULL) {
    // Subsequent calls should just check for accessobj being NULL.
    if (accessobj != NULL) {
      throw AsapError("RGL::SetAtoms called multiple times with accessobj != NULL");
    }
    // SetAtoms should not do anything if called more than once!
    return;
  }

  // The first time SetAtoms is being called some initialization is done.
  if (accessobj != NULL) {
    atoms = accessobj;
    AsapAtoms_INCREF(atoms);
  } else {
    atoms = new NormalAtoms();
  }

  assert(atoms != NULL);
  atoms->Begin(pyatoms);
  nAtoms = atoms->GetNumberOfAtoms();
  nSize = 0;  // Will be set when the NB list is created.  The value 0
              // should create enough trouble to detect if it is used before then.
  atoms->End();

  DEBUGPRINT;
}

// Miscellaneous functions
void RGL::Allocate() {
  USETIMER("RGL::Allocate");
  DEBUGPRINT;
  VERB("Allocate[" << nAtoms << "," << nSize << "]" << flush);

  if (nSize != forces.size() || nAtoms != epot.size()) {
    sigma_p.resize(nAtoms);
    sigma_q.resize(nAtoms);
    epot.resize(nAtoms);
    forces.resize(nSize);
    virials.resize(nSize);
  }

  DEBUGPRINT;
}

long RGL::PrintMemory() const {
  long mem = 0;

  mem += sigma_p.size() * sizeof(double);
  mem += sigma_q.size() * sizeof(double);
  mem += epot.size() * sizeof(double);
  mem += forces.size() * sizeof(Vec);
  mem = (mem + 512*1024) / (512*1024);

  char buffer[500];
  snprintf(buffer, 500,
           "*MEM* RGL %ld MB. [sizeof(int)=%ld sizeof(double)=%ld]",
           mem, (long) sizeof(int), (long) sizeof(double));
  //cerr << buffer << endl;

  if (atoms != NULL)
    mem += atoms->PrintMemory();
  if (neighborlist != NULL)
    mem += neighborlist->PrintMemory();
  return mem;
}

// Neighborlist functions
void RGL::CreateNeighborList() {
  USETIMER("RGL::CreateNeighborList");
  DEBUGPRINT;

  PyAsap_NeighborLocatorObject *nbl;
  nbl = PyAsap_NewNeighborList(atoms, cutoff.end, driftfactor);
  neighborlist = nbl->cobj;
  neighborlist_obj = (PyObject *) nbl;
  neighborlist->UpdateNeighborList();

  DEBUGPRINT;
}

void RGL::UpdateNeighborList()
{
  USETIMER("RGL::UpdateNeighborList");
  DEBUGPRINT;

  if (neighborlist) {
    DEBUGPRINT;
    neighborlist->UpdateNeighborList();
    if ((nAtoms != atoms->GetNumberOfAtoms()) ||
        (nSize - nAtoms != atoms->GetNumberOfGhostAtoms())) {
      nAtoms = atoms->GetNumberOfAtoms();
      nSize = nAtoms + atoms->GetNumberOfGhostAtoms();
      ghostatoms = atoms->HasGhostAtoms();
      Allocate();
    }
  } else {
    DEBUGPRINT;
    VERB("C");
    CreateNeighborList();
    nAtoms = atoms->GetNumberOfAtoms();
    nSize = nAtoms + atoms->GetNumberOfGhostAtoms();
    ghostatoms = atoms->HasGhostAtoms();
    Allocate();
  }
  DEBUGPRINT;
}

void RGL::CheckNeighborList()
{
  USETIMER("RGL::CheckNeighborList");
  DEBUGPRINT;
  VERB(" CheckNBL[");
  assert(atoms != NULL);

  // Update if invalid
  bool update = ((neighborlist == NULL) || neighborlist->IsInvalid());  
  if (!update) {
    DEBUGPRINT;
    VERB("u");
    update = neighborlist->CheckNeighborList();
  }

  //## May communicate
  update = atoms->UpdateBeforeCalculation(update, cutoff.end * (1 + driftfactor));
  if (update) {
    DEBUGPRINT;
    VERB("U");
    UpdateNeighborList();
  }

  DEBUGPRINT;
  VERB("]" << flush);
}

// Get functions
double RGL::GetPotentialEnergy(PyObject *pyatoms) {
  USETIMER("RGL::GetPotentialEnergy");
  DEBUGPRINT;
  VERB(" Energy[");

  double etot = 0.0;
  const vector<double> &energies = GetPotentialEnergies(pyatoms);
  for (int i = 0; i < nAtoms; i++)
    etot += energies[i];
  DEBUGPRINT;
  VERB("]" << endl << flush);
  return etot;
}

const vector<double> &RGL::GetPotentialEnergies(PyObject *pyatoms) {
  USETIMER("RGL::GetPotentialEnergies");
  DEBUGPRINT;
  VERB(" Energies[");

  atoms->Begin(pyatoms);
  CheckNeighborList();
  if (counters.energies != atoms->GetPositionsCounter()) {
    CalculateSigmasAndEnergies();
    counters.energies = atoms->GetPositionsCounter();
  }
  atoms->End();
  DEBUGPRINT;
  VERB("]" << flush);
  return epot;
}

const vector<Vec> &RGL::GetForces(PyObject *pyatoms) {
  USETIMER("RGL::GetForces");
  DEBUGPRINT;
  VERB(" Forces[");

  atoms->Begin(pyatoms);
  CheckNeighborList();
  if (counters.energies != atoms->GetPositionsCounter()) {
    CalculateSigmasAndEnergies();
    counters.energies = atoms->GetPositionsCounter();
  }
  if (counters.forces != atoms->GetPositionsCounter()) {
    CalculateForcesAfterSigmas();
    counters.forces = atoms->GetPositionsCounter();
  }
  atoms->End();
  DEBUGPRINT;
  VERB("]" << endl << flush);
  return forces;
}

const vector<SymTensor> &RGL::GetVirials(PyObject *pyatoms)
{
  USETIMER("RGL::GetVirials");
  DEBUGPRINT;
  VERB(" Virials[");

  atoms->Begin(pyatoms);
  CheckNeighborList();
  if (counters.energies != atoms->GetPositionsCounter()) {
    CalculateSigmasAndEnergies();
    counters.energies = atoms->GetPositionsCounter();
  }
  if (counters.forces != atoms->GetPositionsCounter()) {
    CalculateForcesAfterSigmas();
    counters.forces = atoms->GetPositionsCounter();
  }
  atoms->End();
  DEBUGPRINT;
  VERB("]" << endl << flush);
  return virials;
}


// Calculation required functions
bool RGL::CalcReq_Energy(PyObject *pyatoms) {
  atoms->Begin(pyatoms);
  bool required = (counters.energies != atoms->GetPositionsCounter());
  atoms->End();
  return required;
}

bool RGL::CalcReq_Forces(PyObject *pyatoms) {
  atoms->Begin(pyatoms);
  bool required = (counters.forces != atoms->GetPositionsCounter());
  atoms->End();
  return required;
}

// Calculate functions
void RGL::CalculateSigmasAndEnergies() {
  USETIMER("RGL::CalculateSigmasAndEnergies")
  DEBUGPRINT;
  VERB(" CalcSigEng");

  for (int i = 0; i < nAtoms; i++) {
    sigma_p[i] = 0.0;
    sigma_q[i] = 0.0;
  }

  const asap_z_int *z = atoms->GetAtomicNumbers();
  int maxNeighbors = neighborlist->MaxNeighborListLength();
  vector<int> neighbors(maxNeighbors);
  vector<double> diffs2(maxNeighbors);
  vector<Vec> diffs(maxNeighbors);

  for (int n = 0; n < nAtoms; n++) {
    int i = zmap[z[n]];
    int size = maxNeighbors;
    int nNeighbors = neighborlist->GetNeighbors(n, &neighbors[0], &diffs[0],
                                                &diffs2[0], size);
    for (int k = 0; k < nNeighbors; k++) {
      int m = neighbors[k];
      int j = zmap[z[m]];

      double dr = sqrt(diffs2[k]);

      if (dr < cutoff.start) {
        double dx = dr / param.r0[i][j] - 1.0;
        double dp = param.A[i][j] * exp(-param.p[i][j] * dx);
        double dq = param.qsi2[i][j] * exp(-2.0 * param.q[i][j] * dx);

        sigma_p[n] += dp;
        sigma_q[n] += dq;
        if (m < nAtoms) {
          sigma_p[m] += dp;
          sigma_q[m] += dq;
        }
      } else {
        double dx = dr - cutoff.end;
        double dx3 = dx * dx * dx;
        double dx4 = dx3 * dx;
        double dx5 = dx4 * dx;
        double dp = param.p5[i][j] * dx5 +
                    param.p4[i][j] * dx4 +
                    param.p3[i][j] * dx3;
        double dq = param.q5[i][j] * dx5 +
                    param.q4[i][j] * dx4 +
                    param.q3[i][j] * dx3;

        sigma_p[n] += dp;
        sigma_q[n] += dq * dq;
        if (m < nAtoms) {
          sigma_p[m] += dp;
          sigma_q[m] += dq * dq;
        }
      }
    }

    epot[n] = sigma_p[n] - sqrt(sigma_q[n]);
  }
  DEBUGPRINT;
}

void RGL::CalculateForcesAfterSigmas() {
  USETIMER("RGL::CalculateForcesAfterSigmas")
  DEBUGPRINT;
  VERB(" CalcForAfSig");

  for (int i = 0; i < nSize; i++)
    {
      forces[i][0] = forces[i][1] = forces[i][2] = 0.0;
      virials[i][0] = virials[i][1] = virials[i][2] =
          virials[i][3] = virials[i][4] = virials[i][5] = 0.0;
    }
  const asap_z_int *z = atoms->GetAtomicNumbers();
  int maxNeighbors = neighborlist->MaxNeighborListLength();
  vector<int> neighbors(maxNeighbors);
  vector<double> diffs2(maxNeighbors);
  vector<Vec> diffs(maxNeighbors);

  for (int n = 0; n < nAtoms; n++) {
    int i = zmap[z[n]];
    int size = maxNeighbors;
    int nNeighbors = neighborlist->GetNeighbors(n, &neighbors[0], &diffs[0],
                                                &diffs2[0], size);

    for (int k = 0; k < nNeighbors; k++) {
      int m = neighbors[k];
      int j = zmap[z[m]];

      double dr = sqrt(diffs2[k]);
      double for_p = 0.0;
      double for_q = 0.0;
      if (dr < cutoff.start) {
        double dx = dr / param.r0[i][j] - 1.0;
        for_p = -2.0 * param.A[i][j] * param.p[i][j] *
                exp(-param.p[i][j] * dx) / param.r0[i][j];
        for_q = -param.qsi2[i][j] * param.q[i][j] *
                exp(-2.0 * param.q[i][j] * dx) / param.r0[i][j];
      } else {
        double dx = dr - cutoff.end;
        double dx2 = dx * dx;
        double dx3 = dx2 * dx;
        double dx4 = dx3 * dx;
        double dx5 = dx4 * dx;
        for_p = 2.0 * (5.0 * param.p5[i][j] * dx4 +
                       4.0 * param.p4[i][j] * dx3 +
                       3.0 * param.p3[i][j] * dx2);
        for_q = (param.q5[i][j] * dx5 +
                 param.q4[i][j] * dx4 +
                 param.q3[i][j] * dx3) *
                (5.0 * param.q5[i][j] * dx4 +
                 4.0 * param.q4[i][j] * dx3 +
                 3.0 * param.q3[i][j] * dx2);
      }

      double for_n;
      if (m < nAtoms)
        for_n = (for_p - for_q * (1.0 / sqrt(sigma_q[n]) + 1.0 / sqrt(sigma_q[m]))) / dr;
      else
        for_n = (0.5 * for_p - for_q / sqrt(sigma_q[n])) / dr;
      forces[n] += for_n * diffs[k];
      forces[m] -= for_n * diffs[k];

      for (int alpha = 0; alpha < 3; alpha++) {
        for (int beta = alpha; beta < 3; beta++) {
          double dsx = 0.5 * for_n * diffs[k][alpha] * diffs[k][beta];
          int ii = stresscomp[alpha][beta];
          virials[n][ii] += dsx;
          virials[m][ii] += dsx;
        }
      }
    }
  }

  DEBUGPRINT;
}
