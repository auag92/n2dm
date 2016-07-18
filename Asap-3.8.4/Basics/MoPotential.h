#include "Potential.h"
#include <vector>
#define M 25
using std::vector;

namespace ASAPSPACE {

class Atoms;
class NeighborList2;

class MoPotential : public Potential
{
 public:
  MoPotential(const double *a);
  ~MoPotential();
  
  virtual string GetName() const {return "MoPotential";}
  
  void SetAtoms(Atoms *a); 
  //virtual void SetPar(const double *a);
  double GetPotentialEnergy();
  const Vec *GetCartesianForces();
  const symTensor *GetStresses(const Vec *momenta = 0);
  void GetStress(double stress[6], const Vec *momenta = 0);
  const double *GetPotentialEnergies();
  void CheckNeighborLists();
  void CalculateRhos();

  double GetCutoffRadius() const {return rCut;}
  double GetLatticeConstant() const {return latticeConstant;}
  int GetNumberOfAtoms() const {return nAtoms;}
  void UpdateSuperCell(const SuperCell *newSuperCell) {}

  /// This potential can be used in parallel simulations
  virtual bool Parallelizable() const {return true;}

private:
  void Allocate();
  double CalculateEnergyAndEnergies();
  void CalculateForcesAndStresses(Vec forces[]);
  //Atoms *atoms;
  NeighborList2 *neighborList;
  int nAtoms;
  int nSize;
  
  double Expm(double a, double x);
  double dExpm(double a, double x);
  double func(double x);
  double dfunc(double x);
  void funcPointer(double x, double *value, double *dvalue);
  double Density(double r);
  double dDensity(double r);
  void DensityPointer(double r, double *dens, double *ddens);
  double PairPotential(double r);
  double dPairPotential(double r);
  void PairPotentialPointer(double r, double *pair, double *dpair);
  double PairPotentialRef(double a);
  double dPairPotentialRef(double a);
  double EmbeddingEnergy(double a);
  double dEmbeddingEnergy(double a);
  double LatticeConstant(double rho);
  double dLatticeConstant(double rho);

  double rCut, latticeConstant, alpha, c[2], norm, coef[M], funccoef[13], latticefunclist[8], rhoMin, rhoMinMin, rhoMax, kmin, kmax, K, expArCut, E0, dE0, gamma;
  double rInCut,dpair0,pair0;
  int alphaMax, coefMax, funccoefMax;
  vector<double> atomicEnergies;
  vector<Vec> potforces;
  vector<Vec> forces; // GetCartesianForces returns a pointer into this one.
  vector<double> stresses;  // GetStresses returns pointer into this one.
  vector<double> potenergies;
  vector<double> potstresses;
  vector<double> lat;
  vector<double> dens0;
  vector<double> Factor;
  vector<double> rho1;
/*   vector<double> xdens; */
/*   vector<double> ydens; */
/*   vector<double> zdens; */
  vector<double> Xdens;  
  vector<double> rho2;
  vector<double> dens2;
/*   vector<double> xxdens; */
/*   vector<double> yydens; */
/*   vector<double> zzdens; */
/*   vector<double> xydens; */
/*   vector<double> yzdens; */
/*   vector<double> zxdens;   */
  vector<double> XXdens;  

  double unnormalizedstress[6];  // The total stress before division by volume
  double totalvolume;            // Total volume during stresscalculation

// A structure of counters to check if recalculations are necessary.
  struct {
    int nblist;
    int energies;
    int forces;
    int stresses;
  } counters;
};

} // end namespace
