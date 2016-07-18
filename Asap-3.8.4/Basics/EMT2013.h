/*
 * EMT2013.h, An update to EMT.h introducing changes to three methods,
 * sigma_batch, force_batch and CalculateEnergiesAfterSigmas
 *
 *  Created on: Apr 29, 2011
 *      Author: s072162
 */

#ifndef _EMT2013_H
#define _EMT2013_H

//// Uses everything from EMT.h except the three changes
#include "EMT.h"

namespace ASAPSPACE {

class EMT2013 : public EMT
{
public:
  //// Creates an EMT2013 potential
  EMT2013(PyObject *self, PyObject *params, bool no_new_elements);

  //// Destructor.
  virtual ~EMT2013();

  /// Return a copy (new reference) of the parameter dictionary
  PyObject *GetParameterDict() const;

  /* Functions changed from the earlier EMT-version. */
protected:
  // Parts of the Energy calculation
  virtual void CalculateEnergiesAfterSigmas(bool calcEpot);

  virtual void sigma_batch(int *self, int *other, Vec rnb[],
      double *sq_dist, int zs, int zo, int n, bool calculatesigma2,
      bool partialupdate = false);

  // Part of the force calculation
  virtual void force_batch(const int *self, const int *other, const Vec rnb[],
      const double sq_dist[], const double dEdss[],
      const double dEdso[], int zs, int zo, int n);

  virtual void CreateNeighborList();
  virtual void CalculateIDs();

  void InitParameters();

  void GetListOfElements(set<int> &elements);
  emt_parameters *ExtractParameters(int z);
  double ExtractParam_double(PyObject *params_for_z, const char *param_name);
  void CalculateCutoffDistances();
  void CalculateChi();

protected:
  // New parameters
  TinyMatrix<double> rcut2;
  TinyMatrix<double> rcut2_NB;
  PyObject *param_obj;
  static const float beta;
  bool no_new_elements;
  int id_assigned;  // Used in CalculateIDs
};

} // end namespace

#endif /* _EMT2013_H */
