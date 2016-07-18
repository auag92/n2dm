#ifndef _RGL_H
#define _RGL_H

#include "AsapPython.h"
#include "Asap.h"
#include "Vec.h"
#include "SymTensor.h"
#include "Potential.h"
#include "TinyMatrix.h"
#include <vector>
using std::vector;

namespace ASAPSPACE {

class Atoms;

class RGL : public Potential {
  public:
    RGL(PyObject *self,
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
        const double &rce);

    ~RGL();

    void SetAtoms(PyObject *pyatoms, Atoms* accessobj = NULL);

    // Get functions
    const vector<double> &GetPotentialEnergies(PyObject *pyatoms);
    double GetPotentialEnergy(PyObject *pyatoms);
    const vector<Vec> &GetForces(PyObject *pyatoms);
    const vector<SymTensor> &GetVirials(PyObject *pyatoms);

    string GetName() const {return "RGL";}
    PyObject *GetNeighborList() const {return neighborlist_obj;}
    double GetCutoffRadius() const {return cutoff.end;}
    int GetNumberOfAtoms() const {return nAtoms;}
    double GetLatticeConstant() const {return 0.0;}

    // Calculation required functions
    bool CalcReq_Energy(PyObject *pyatoms);
    bool CalcReq_Forces(PyObject *pyatoms);
    bool CalcReq_Virials(PyObject *pyatoms) {return CalcReq_Forces(pyatoms);}

    // Miscellaneous functions
    bool Parallelizable() const {return true;}
    long PrintMemory() const;

  protected:
    void Allocate();

    // Neighborlist functions
    void CreateNeighborList();
    void UpdateNeighborList();
    void CheckNeighborList();

    // Calculate functions
    void CalculateSigmasAndEnergies();
    void CalculateForcesAfterSigmas();
    void CalculateStressesAfterForces(bool virialonly);

    // Data structures and variables
    // Atoms *atoms;

    int nElements;
    int nAtoms;
    int nSize;
    double driftfactor;
    vector<int> zmap;

    NeighborLocator *neighborlist;
    PyObject *neighborlist_obj;
    bool ghostatoms;

    struct {
      int energies;
      int forces;
    } counters;

    struct {
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
    } param;

    struct {
      double start;
      double end;
    } cutoff;

    vector<double> sigma_p;
    vector<double> sigma_q;
    vector<double> epot;
    vector<Vec> forces;
    vector<SymTensor> virials;
};

} // end namespace

#endif // _RGL_H

