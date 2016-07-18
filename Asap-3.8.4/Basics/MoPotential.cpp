#include "MoPotential.h"
#include "NeighborList2.h"
#include "Atoms.h"
#include "Vec.h"
#include "GhostAtoms.h"
#include "GhostPotential.h"
// #include "Exception.h"
#include "SuperCell.h"
// #include "ParallelAtoms.h"
#include <math.h>
#include <assert.h>
#include <stdio.h>
#define PI 3.1415926543
#define M 25

using std::cerr;
using std::endl;
using std::flush;

#if 0
#define DEBUGPRINT cerr << "Now in " << __FILE__ << " line " << __LINE__ << endl << flush;
#else
#define DEBUGPRINT
#endif

#ifdef SUN_SQRT_LITERAL_HACK
// SunOS compiler: sqrt does not work with a literal number (cannot decide if
// it is a double or a long double).
#define sqrt(x) sqrt((double) x)
#endif

const static int stresscomp[3][3] = {{0, 5, 4}, {5, 1, 3}, {4, 3, 2}};
//const static int stresscomp[3][3] = {{0, 1, 2}, {1, 3, 4}, {2, 4, 5}};

MoPotential::MoPotential(const double *a)
{
  DEBUGPRINT
  rCut = 4.85;
  double Rin = 1.85;
  K = PI/(rCut-Rin);
  latticeConstant = 3.1872;
  alpha = a[0];  
  expArCut = exp(-alpha*rCut);
  gamma = 0.91;
  E0 = exp(-1.561019*latticeConstant*gamma)*(5632.825549-2211.468843*latticeConstant*gamma);
  dE0 = exp(-1.561019*latticeConstant*gamma)*(-11004.418173+3452.145520*latticeConstant*gamma);
  c[0] = a[1];
  c[1] = a[2];
  coefMax = 0;
  for(int i = 0;i < M;i++)
    {
    coef[i] = a[i+3];
    if (coef[i]!=0.0) coefMax = i+1;
    }
  funccoef[0] = 1.0;
  funccoefMax = 1;
  for(int i = 1;i < 13;i++)
    {
    funccoef[i] = a[i+2+M];
    if (funccoef[i]!=0.0) funccoefMax = i+1;
    }
  for(int i = 0;i < 8;i++)
    {
    latticefunclist[i] = a[i+15+M];
    }
  norm = 8*Density(latticeConstant*sqrt(3.0)/2.0)+6*Density(latticeConstant)+12*Density(latticeConstant*sqrt(2.0));
  rhoMin = (8*Density(3.4*sqrt(3.0)/2.0)+6*Density(3.4)+12*Density(3.4*sqrt(2.0)))/norm;
  rhoMinMin = (8*Density(4.1)+6*Density(4.1*2/sqrt(3))+12*Density(4.1*sqrt(2.0)))/norm;
  rhoMax = (8*Density(2.9*sqrt(3.0)/2.0)+6*Density(2.9)+12*Density(2.9*sqrt(2.0)))/norm;
  kmin = exp(-3.4/rhoMin/dLatticeConstant(rhoMin))*norm*rhoMin;
  kmax = exp(-2.9/rhoMax/dLatticeConstant(rhoMax))*norm*rhoMax; 
  //cerr << kmin << " " << kmax << endl;
  neighborList = 0;
  //cerr << rhoMin << " " << rhoMax << endl;
  //cerr << LatticeConstant(rhoMin+.001) << " " << LatticeConstant(1.) << " " << LatticeConstant(rhoMax-.001) << endl;  
  counters.nblist = counters.energies = counters.forces = counters.stresses = 0;      
  double *pt_e, *pt_df;
  pt_e = &pair0;
  pt_df = &dpair0;
  rInCut = 0.1;//2.1
  PairPotentialPointer(2.2, pt_e, pt_df);//2.2
  rInCut = 2.2;//2.2
  double eold,dfold;
  pt_e = &eold;
  pt_df = &dfold;
  PairPotentialPointer(1.99, pt_e, pt_df);
  for (int i = 400;i<400;i++)
    {
      double e,df;
      pt_e = &e;
      pt_df = &df;
      PairPotentialPointer(i/100., pt_e, pt_df);
      cerr << i/100. << " " << e << " " << df << " " << (e-eold)/0.01 << " " << (df+dfold)/2 << endl;
      eold = e;
      dfold = df;
    }
  for (int i = 34950;i<34050;i++)
    {
      double rho = (8*Density(i/10000.*sqrt(3.0)/2.0)+6*Density(i/10000.)+12*Density(i/10000.*sqrt(2.0)))/norm;
      double lat = LatticeConstant(rho);
      cerr << i/10000. << " " << lat << " " << dLatticeConstant(rho) << " " << dDensity(i/20000.*sqrt(3))/norm*(dEmbeddingEnergy(lat)-dPairPotentialRef(lat))*dLatticeConstant(rho) << endl;
    }
  for (int i = 29950;i<29050;i++)
    {
      double rho = (8*Density(i/10000.*sqrt(3.0)/2.0)+6*Density(i/10000.)+12*Density(i/10000.*sqrt(2.0)))/norm;
      double lat = LatticeConstant(rho);
      cerr << i/10000. << " " << lat << " " << dLatticeConstant(rho) << " " << dDensity(i/20000.*sqrt(3))/norm*(dEmbeddingEnergy(lat)-dPairPotentialRef(lat))*dLatticeConstant(rho) << endl;
    }
}

MoPotential::~MoPotential()
{
  DEBUGPRINT;
  if (neighborList) 
    {
      delete neighborList;
    }
}

void MoPotential::SetAtoms(Atoms *a)
{
  DEBUGPRINT;
  atoms = a;
  DEBUGPRINT;
}

void MoPotential::Allocate()
{
  DEBUGPRINT;
  // In case this is a QuasiContinuum calculations, where "atoms" points
  // to a "QCAtoms" object, we want "nAtoms" to be the number of atomic
  // atoms and not the number of atoms plus the number of nodes. That is
  // why we ask the "Atoms" base class for the number of atoms:
  nAtoms = atoms->GetNumberOfRealAtoms();
  nSize = nAtoms;
  // Do we have any ghosts?
  GhostAtoms *ghostAtoms = dynamic_cast<GhostAtoms *>(atoms);
  if (ghostAtoms != 0)
    nSize += ghostAtoms->GetNumberOfGhosts();    // Yes!
  atomicEnergies.resize(nSize);
  potforces.resize(nSize);
  forces.resize(nSize);
  stresses.resize(6*nSize);
  potstresses.resize(6*nSize);
  lat.resize(nSize);
  dens0.resize(nSize);
  Factor.resize(nSize);
  rho1.resize(nSize);
//   xdens.resize(nSize);
//   ydens.resize(nSize);
//   zdens.resize(nSize);
  Xdens.resize(3*nSize);
  rho2.resize(nSize);
  dens2.resize(nSize);
//   xxdens.resize(nSize);
//   yydens.resize(nSize);
//   zzdens.resize(nSize);
//   xydens.resize(nSize);
//   yzdens.resize(nSize);
//   zxdens.resize(nSize);
  XXdens.resize(6*nSize);
}

void MoPotential::CheckNeighborLists()
{
    DEBUGPRINT;
    if (counters.nblist == atoms->GetChangeCounter())
      return;
    counters.nblist = atoms->GetChangeCounter();
    if (neighborList)
      {   
	DEBUGPRINT;    
        bool reallocate = neighborList->CheckAndUpdateNeighborList();
	DEBUGPRINT;
        if (reallocate)
            Allocate();
      }
    else
      {
	DEBUGPRINT;
        neighborList = new NeighborList2(atoms, rCut);
	DEBUGPRINT;
        neighborList->CheckAndUpdateNeighborList();
	DEBUGPRINT;
        Allocate();
      }
    DEBUGPRINT;
}

double MoPotential::GetPotentialEnergy()
{
  DEBUGPRINT;
  CheckNeighborLists();
  DEBUGPRINT;
  return CalculateEnergyAndEnergies();
}

const double *MoPotential::GetPotentialEnergies()
{
  DEBUGPRINT;
  CheckNeighborLists();
  CalculateEnergyAndEnergies();
  DEBUGPRINT;
  return &(atomicEnergies[0]);
}

const Vec *MoPotential::GetCartesianForces()
{
  DEBUGPRINT;
  CheckNeighborLists();
  memset(&(forces[0]), 0, nAtoms * sizeof(Vec));
  CalculateForcesAndStresses(&(forces[0]));
  DEBUGPRINT;
  return &(forces[0]);
}

const symTensor *MoPotential::GetStresses(const Vec *momenta)
{
  //  if (momenta)
  //   cerr << "WARNING: Dynamic part of stress is not calculated!" << endl;
  CheckNeighborLists();
  symTensor *stresses = (symTensor *) &(this->stresses[0]);
  Vec *dummyForces = new Vec[nAtoms];
  memset(stresses, 0, nAtoms * sizeof(symTensor));
  if (counters.stresses != atoms->GetChangeCounter())
    {
      CalculateForcesAndStresses(dummyForces);
    }
  double (*potstresses)[6];
  potstresses = (double (*)[6]) &(this->potstresses[0]);
  for (int a1 = 0; a1 < nAtoms; a1++)
    for(int i = 0;i<6;i++)
      stresses[a1][i] = potstresses[a1][i];
  delete [] dummyForces;
  // from EMT.cpp
  //--------------------------------------------------------
  for (int i = 0; i < 6; i++)
    unnormalizedstress[i] = 0.0;
  totalvolume = 0.0;
  double vol = 3.1877 * 3.1877 * 3.1877 * 0.5;
  double invvol = 1.0 / vol;
  totalvolume = nAtoms*vol;
  for (int i = 0; i < nAtoms; i++)
    for (int alpha = 0; alpha < 3; alpha++)
      for (int beta = alpha; beta < 3; beta++)
	{
	  int j = stresscomp[alpha][beta];
	  if (momenta)
	    stresses[i][j] -= momenta[i][alpha] * momenta[i][beta] / 95.94;
	  unnormalizedstress[j] += stresses[i][j];
	  stresses[i][j] *= invvol;
	}
  //--------------------------------------------------------  
  return stresses;
}

// void MoPotential::CalculateStress(double stress[6], Vec *momenta)
// {
//   if (momenta)
//     cerr << "WARNING: Dynamic part of stress is not calculated!" << endl;
//   CheckNeighborLists();
//   Vec *dummyForces = new Vec[nAtoms];
//   double (*stresses)[6] = new double[nAtoms][6];
//   memset(stresses, 0, nAtoms * sizeof(symTensor));
//   CalculateForcesAndStresses(dummyForces, stresses);
//   for (int k = 0; k < 6; k++)
//     {
//       stress[k] = 0;
//       for (int a = 0; a < nAtoms; a++)
//         stress[k] += stresses[a][k];
//     }
//   delete [] stresses;
//   delete [] dummyForces;
// }

void MoPotential::GetStress(double stress[6], const Vec *momenta)
{
  double vol;
  DEBUGPRINT;

  // Calculate the atomic stresses. Their sum is stored in unnormalizedstress.
  GetStresses(momenta);
  
  // Always use the volume of the supercell to normalize the total stress.
  // The alternative, using the sum of atomic volumes when the supercell
  // volume is ill-defined due to open boundary conditions, will fail for
  // parallel simulations.
  vol = atoms->GetSuperCell()->GetVolume();
    
  // Normalize the total stress with the volume.
  for (int i = 0; i < 6; i++)
    stress[i] = (stress[i] + unnormalizedstress[i]) / vol;
  //stress[i] = (stress[i] + unnormalizedstress[i]) / totalvolume;
  // WHY ADD TO OLD VALUE ????
  DEBUGPRINT;
}

double MoPotential::CalculateEnergyAndEnergies()
{
  if (counters.energies != atoms->GetChangeCounter())
    {
      counters.energies = atoms->GetChangeCounter();
      CalculateRhos();
    }
  double energy = 0.0;
  for (int a = 0; a < nAtoms; a++)
    energy +=  atomicEnergies[a];
  DEBUGPRINT;
  return energy;
}

double MoPotential::Expm(double a, double x)
{
  if (x<rCut)
    return exp(-a*x)-(a*(rCut-x)+1)*expArCut;
//return exp(-a*x)-(a*(rCut-x)+1)*exp(-a*rCut);
  else
    return 0.;
}

double MoPotential::dExpm(double a, double x)
{
  if (x<rCut)
    return a*(expArCut-exp(-a*x));
//return a*(exp(-a*rCut)-exp(-a*x));
  else
    return 0.;
}

double MoPotential::func(double x)
{
  double funcvalue = 0.0;
  if (x<rCut)
    for (int i = 0;i < funccoefMax;i++)
      {
	funcvalue += funccoef[i]/(K*(i+1))*(1-cos((K*(i+1))*(rCut-x)));
      }
  return funcvalue;
}

double MoPotential::dfunc(double x)
{
  double dfuncvalue = 0.0;
  if (x<rCut)  
    for (int i = 0;i < funccoefMax;i++)
      {
	dfuncvalue += -funccoef[i]*sin((K*(i+1))*(rCut-x));
      }
  return dfuncvalue;
}

void MoPotential::funcPointer(double r, double *value, double *dvalue)
{
  double cosntheta,sinntheta,costheta,sintheta,tmpcos=1,tmpsin=0;
  *value = 0.0;
  *dvalue = 0.0;
  if (r<rCut)
    {
      costheta = cos(K*(rCut-r));
      sintheta = sin(K*(rCut-r));
      for (int i = 0;i < funccoefMax;i++)
	{
	  cosntheta = tmpcos*costheta-tmpsin*sintheta;
	  sinntheta = tmpsin*costheta+tmpcos*sintheta;  
	  *value += funccoef[i]/(K*(i+1))*(1-cosntheta);
	  *dvalue += -funccoef[i]*sinntheta;
	  tmpcos = cosntheta;
	  tmpsin = sinntheta;
      }
    }
}

double MoPotential::Density(double r)
{
  return Expm(alpha,r);
}

double MoPotential::dDensity(double r)
{
  return dExpm(alpha,r);
}

void MoPotential::DensityPointer(double x, double *dens, double *ddens)
{
  if (x<rCut)
    {
      double expX = exp(-alpha*x);
      *dens = expX-(alpha*(rCut-x)+1)*expArCut;
      *ddens = alpha*(expArCut-expX);
    }
  else
    {
      *dens = 0.0;
      *ddens = 0.0;
    }
}

double MoPotential::PairPotential(double r)
{
  double pair = 0.0;
  double cosntheta,sinntheta,costheta,sintheta,tmpcos=1,tmpsin=0;
  if (r<rCut)
    {
      if (r>rInCut)
	{
	  costheta = cos(K*(rCut-r));
	  sintheta = sin(K*(rCut-r));
	  for (int i = 0;i < coefMax;i++)
	    {
	      cosntheta = tmpcos*costheta-tmpsin*sintheta;
	      sinntheta = tmpsin*costheta+tmpcos*sintheta;
	      //pair += coef[i]/(K*(i+1))*(1-cos((K*(i+1))*(rCut-r)));
	      pair += coef[i]/(K*(i+1))*(1-cosntheta);
	      tmpcos = cosntheta;
	      tmpsin = sinntheta;
	    }
	}
      else
	{
	  //	  cerr << r << " <  " << rInCut << endl;
	  return (r-rInCut)*dpair0 + pair0;
	}
    }
  return pair;
}

double MoPotential::dPairPotential(double r)
{
  double dpair = 0.0;
  double cosntheta,sinntheta,costheta,sintheta,tmpcos=1,tmpsin=0;
  if (r<rCut)
    {
      if (r>rInCut)
	{
	  costheta = cos(K*(rCut-r));
	  sintheta = sin(K*(rCut-r));
	  for (int i = 0;i < coefMax;i++)
	    {
	      cosntheta = tmpcos*costheta-tmpsin*sintheta;
	      sinntheta = tmpsin*costheta+tmpcos*sintheta;
	      //dpair += -coef[i]*sin((K*(i+1))*(rCut-r));
	      dpair += -coef[i]*sinntheta;
	      tmpcos = cosntheta;
	      tmpsin = sinntheta;
	    }
	}
      else
	{
	  //	  cerr << r << " <  " << rInCut << endl;
	  return dpair0;
	}
    }
  return dpair;
}

void MoPotential::PairPotentialPointer(double r, double *pair, double *dpair)
{
  *pair = 0.0;
  *dpair = 0.0;
  double cosntheta,sinntheta,costheta,sintheta,tmpcos=1,tmpsin=0;
  if (r<rCut)
    {
      if (r>rInCut)
	{
	  costheta = cos(K*(rCut-r));
	  sintheta = sin(K*(rCut-r));
	  for (int i = 0;i < coefMax;i++)
	    {
	      cosntheta = tmpcos*costheta-tmpsin*sintheta;
	      sinntheta = tmpsin*costheta+tmpcos*sintheta;
	      //pair += coef[i]/(K*(i+1))*(1-cos((K*(i+1))*(rCut-r)));
	      *pair += coef[i]/(K*(i+1))*(1-cosntheta);
	      *dpair += -coef[i]*sinntheta;
	      tmpcos = cosntheta;
	      tmpsin = sinntheta;
	    }
	}
      else
	{
	  //	  cerr << r << " <  " << rInCut << endl;
	  *pair = (r-rInCut)*dpair0 + pair0;
	  *dpair = dpair0;	  
	}
    }
}

double MoPotential::PairPotentialRef(double a)
{
  return 8*PairPotential(a*sqrt(3.0)/2.0)+6*PairPotential(a)+12*PairPotential(a*sqrt(2.0));
}

double MoPotential::dPairPotentialRef(double a)
{
  return 4*sqrt(3.0)*dPairPotential(a*sqrt(3.0)/2.0)+6*dPairPotential(a)+12*sqrt(2.0)*dPairPotential(a*sqrt(2.0));
}

double MoPotential::EmbeddingEnergy(double a)
{
  //return exp(-1.848042798*a)*(12641.649765-4845.9460549*a);
  //return exp(-1.561019*a)*(5632.825549-2211.468843*a);
  if (a>latticeConstant*gamma)
    return exp(-1.561019*a)*(5632.825549-2211.468843*a);
  else
    return E0+(a-latticeConstant*gamma)*dE0+32*(a-latticeConstant*gamma)*(a-latticeConstant*gamma);
}

double MoPotential::dEmbeddingEnergy(double a)
{
  //return exp(-1.848042798*a)*(-28208.255859+8955.5157065*a);
  //return exp(-1.561019*a)*(-11004.418173+3452.145520*a);
  if (a>latticeConstant*gamma)
    return exp(-1.561019*a)*(-11004.418173+3452.145520*a);
  else
    return dE0+64*(a-latticeConstant*gamma); 
}

double MoPotential::LatticeConstant(double rho)
{
  if (rho >= rhoMin && rho <= rhoMax)
    {
      //return latticefunclist[0] + (latticefunclist[1]+latticefunclist[2]*rho+latticefunclist[3]*rho*rho+latticefunclist[4]*rho*rho*rho)*log(rho) + rho*latticefunclist[5] + 1/rho*latticefunclist[6] + 1/rho/rho*latticefunclist[7];
      return latticefunclist[0] + (latticefunclist[1]+rho*(latticefunclist[2]+rho*(latticefunclist[3]+rho*latticefunclist[4])))*log(rho) + rho*latticefunclist[5] + 1/rho*(latticefunclist[6]+1/rho*latticefunclist[7]);
    }
  else if (rho > rhoMax)
    {
      return 2.9*log(norm*rho/kmax)/log(norm*rhoMax/kmax);
    } 
  else if (rho < rhoMin && rho >= rhoMinMin)
    {
      return 3.4*log(norm*rho/kmin)/log(norm*rhoMin/kmin);
    }
  else 
    {
      return 3.4*log(norm*rhoMinMin/kmin)/log(norm*rhoMin/kmin);
    }
}

double MoPotential::dLatticeConstant(double rho)
{
  if (rho >= rhoMin && rho <= rhoMax)
    {
      //return latticefunclist[1]/rho+latticefunclist[2]*(log(rho)+1)+latticefunclist[3]*(2*rho*log(rho)+rho)+latticefunclist[4]*(3*rho*rho*log(rho)+rho*rho) + latticefunclist[5] - 1/rho/rho*latticefunclist[6] - 2/rho/rho/rho*latticefunclist[7];
      return (latticefunclist[2]+rho*(2*latticefunclist[3]+3*rho*latticefunclist[4]))*log(rho) + latticefunclist[2] + (latticefunclist[3]+ rho*latticefunclist[4])*rho + latticefunclist[5] + 1/rho*(latticefunclist[1] -1/rho*(latticefunclist[6]+2/rho*latticefunclist[7]));
    }
  else if (rho > rhoMax)
    {
      return 2.9/log(norm*rhoMax/kmax)/rho;
    }
  else if (rho < rhoMin && rho >= rhoMinMin)
    {
      return 3.4/log(norm*rhoMin/kmin)/rho;
    }
  else  
    {
      return 0.0;
    } 
}

void MoPotential::CalculateForcesAndStresses(Vec *forces)
{
  //cerr << "void MoPotential::CalculateForcesAndStresses()" << endl;
  // we assume that the forces have been zeroed.
  // They are if they come from PyArray_FromDims.
  if (counters.energies != atoms->GetChangeCounter())// || stress)
    {
      counters.energies = atoms->GetChangeCounter();
      CalculateRhos();
    }
  const int maxSize = 100;
  int neighbors[maxSize];
  Vec diffs[maxSize];
  double diffs2[maxSize];
  double r[maxSize];
  double invr[maxSize];
  double X[maxSize][3];
  Vec drdR[maxSize];
  double dX[maxSize][3][3];
  double dFX[maxSize][3];
  double tmpFactor[maxSize];
  double (*potstresses)[6];
  potstresses = (double (*)[6]) &(this->potstresses[0]);
  double (*Xdens)[3];
  Xdens  = (double (*)[3]) &(this->Xdens[0]);
  double tmpXdens[maxSize][3];
  double (*XXdens)[6];
  XXdens  = (double (*)[6]) &(this->XXdens[0]);
  double tmpXXdens[maxSize][6];
  double tmpdens2[maxSize];
  double dens[maxSize];
  double ddens[maxSize];
  double tmp[maxSize];
  double drdRddensX[maxSize][3];
  double dXdens[maxSize][3][3];
  DEBUGPRINT;
  for (int a1 = 0; a1 < nAtoms; a1++)
    {
      forces[a1] += potforces[a1];
      int size = maxSize;
      int n = neighborList->
	GetNeighbors(a1, neighbors, diffs, diffs2, size);
      assert(size >= 0);
      for (int m = 0; m < n; m++)
	{
	  r[m] = sqrt(diffs2[m]);
	  invr[m] = 1/r[m];
	}
      for (int m = 0; m < n; m++)
	{
	  drdR[m]  = -diffs[m]*invr[m];
	  //x[m] = -drdR[m][0];//diffs[m][0]*invr[m];
	  //y[m] = -drdR[m][1];//diffs[m][1]*invr[m];
	  //z[m] = -drdR[m][2];//diffs[m][2]*invr[m];
	  for (int i = 0; i < 3; i++)
	    X[m][i] = -drdR[m][i];
	  //dx[m] = Vec(-1.,0.,0.)*invr[m] - drdR[m]*x[m]*invr[m];//+ diffs[m]*x[m]*invr[m]*invr[m];
	  //dy[m] = Vec(0.,-1.,0.)*invr[m] - drdR[m]*y[m]*invr[m];// + diffs[m]*y[m]*invr[m]*invr[m];
	  //dz[m] = Vec(0.,0.,-1.)*invr[m] - drdR[m]*z[m]*invr[m];// + diffs[m]*z[m]*invr[m]*invr[m];
	  for (int i = 0; i < 3; i++)
	    for (int j = i; j < 3; j++)
	      {
		dX[m][i][j] = X[m][i]*X[m][j]*invr[m];
		dX[m][j][i] = dX[m][i][j]; 
	      }
	  for (int i = 0; i < 3; i++)
	    dX[m][i][i] -= invr[m];
	}
      for (int m = 0; m < n; m++)
	{
	  int a2 = neighbors[m];
	  tmpFactor[m] = Factor[a1] + Factor[a2];
	  //tmpxdens[m] = xdens[a1] - xdens[a2];
	  //tmpydens[m] = ydens[a1] - ydens[a2];
	  //tmpzdens[m] = zdens[a1] - zdens[a2];
	  for (int i = 0; i < 3; i++)
	    tmpXdens[m][i] = Xdens[a1][i] - Xdens[a2][i];
	  //tmpxxdens[m] = xxdens[a1] + xxdens[a2];
	  //tmpyydens[m] = yydens[a1] + yydens[a2];
	  //tmpzzdens[m] = zzdens[a1] + zzdens[a2];
	  //tmpxydens[m] = xydens[a1] + xydens[a2];
	  //tmpyzdens[m] = yzdens[a1] + yzdens[a2];
	  //tmpzxdens[m] = zxdens[a1] + zxdens[a2];
	  for (int i = 0; i < 6; i++)
	    tmpXXdens[m][i] = XXdens[a1][i] + XXdens[a2][i];
	  tmpdens2[m] = dens2[a1] + dens2[a2];
	}
      for (int m = 0; m < n; m++)
	{
	  ddens[m] = dDensity(r[m])/norm;	
	}
      for (int m = 0; m < n; m++)
	for (int i = 0; i < 3; i++)
	  dFX[m][i] = X[m][i]*ddens[m]*tmpFactor[m];
//       for (int m = 0; m < n; m++)
// 	{
// 	  dF[m] = -drdR[m]*ddens[m]*tmpFactor[m];
// 	}
      //for (int m = 0; m < n; m++)
      //	cerr << dF[m][0] - dFX[m][0] << " " << dF[m][1] - dFX[m][1] << " " << dF[m][2] - dFX[m][2] << endl; 
      for (int m = 0; m < n; m++)
	{
	  double *pt_dens, *pt_ddens;
	  pt_dens = &dens[m];
	  pt_ddens = &ddens[m];
	  funcPointer(r[m], pt_dens, pt_ddens);
	  dens[m] = *pt_dens*c[0];
	  ddens[m] = *pt_ddens*c[0];
	}
      for (int m = 0; m < n; m++)
	{
	  for (int i = 0; i < 3; i++)
	    for (int j = 0; j < 3; j++)
	      dFX[m][i] -= (dX[m][j][i]*dens[m] - X[m][i]*X[m][j]*ddens[m])*tmpXdens[m][j]*2; 
	  //dF[m] += -(dx[m]*dens[m] + drdR[m]*x[m]*ddens[m])*tmpxdens[m]*2 - (dy[m]*dens[m] + drdR[m]*y[m]*ddens[m])*tmpydens[m]*2 - (dz[m]*dens[m] + drdR[m]*z[m]*ddens[m])*tmpzdens[m]*2;
	}
      for (int m = 0; m < n; m++)
	{
	  tmp[m] = dens[m]*c[1]/c[0];
	}
      for (int m = 0; m < n; m++)
	{
	  dens[m] = tmp[m];	
	}
      for (int m = 0; m < n; m++)
	{
	  tmp[m] = ddens[m]*c[1]/c[0];
	}
      for (int m = 0; m < n; m++)
	{
	  ddens[m] = tmp[m];	
	}
      for (int m = 0; m < n; m++)
	for (int i = 0; i < 3; i++)
	  {
	    drdRddensX[m][i] = -X[m][i]*ddens[m];
	    for (int j = 0; j < 3; j++)
	      dXdens[m][i][j] = dX[m][i][j]*dens[m];
	  }
      for (int m = 0; m < n; m++)
	{
	  //Vec drdRddens = drdR[m]*ddens[m];
	  //Vec dxdens = dx[m]*dens[m];
	  //Vec dydens = dy[m]*dens[m];
	  //Vec dzdens = dz[m]*dens[m];
	  //dF[m] += -(dxdens*x[m]*2 + drdRddens*x[m]*x[m])*tmpxxdens[m]*2 - (dydens*y[m]*2 + drdRddens*y[m]*y[m])*tmpyydens[m]*2 - (dzdens*z[m]*2 + drdRddens*z[m]*z[m])*tmpzzdens[m]*2 - (dxdens*y[m] + dydens*x[m] + drdRddens*x[m]*y[m])*tmpxydens[m]*4 - (dydens*z[m] + dzdens*y[m] + drdRddens*y[m]*z[m])*tmpyzdens[m]*4 - (dzdens*x[m] + dxdens*z[m] + drdRddens*z[m]*x[m])*tmpzxdens[m]*4 + (drdRddens)*tmpdens2[m]*(2/3.); 
	  for (int i = 0; i < 3; i++)
	    {
	      for (int j = 0; j < 3; j++)
		{
		  int k = (j+1)%3;
		  dFX[m][i] += -(dXdens[m][j][i]*X[m][j]*2 + drdRddensX[m][i]*X[m][j]*X[m][j])*tmpXXdens[m][j]*2  - (dXdens[m][j][i]*X[m][k] + dXdens[m][k][i]*X[m][j] + drdRddensX[m][i]*X[m][j]*X[m][k])*tmpXXdens[m][j+3]*4;
		  //dFX[m][i] += -(dXdens[m][0][i]*X[m][0]*2 + drdRddensX[m][i]*X[m][0]*X[m][0])*tmpXXdens[m][0]*2 - (dXdens[m][1][i]*X[m][1]*2 + drdRddensX[m][i]*X[m][1]*X[m][1])*tmpXXdens[m][1]*2 - (dXdens[m][2][i]*X[m][2]*2 + drdRddensX[m][i]*X[m][2]*X[m][2])*tmpXXdens[m][2]*2 - (dXdens[m][0][i]*X[m][1] + dXdens[m][1][i]*X[m][0] + drdRddensX[m][i]*X[m][0]*X[m][1])*tmpXXdens[m][3]*4 - (dXdens[m][1][i]*X[m][2] + dXdens[m][2][i]*X[m][1] + drdRddensX[m][i]*X[m][1]*X[m][2])*tmpXXdens[m][4]*4 - (dXdens[m][2][i]*X[m][0] + dXdens[m][0][i]*X[m][2] + drdRddensX[m][i]*X[m][2]*X[m][0])*tmpXXdens[m][5]*4;
		}
	      dFX[m][i] +=(drdRddensX[m][i])*tmpdens2[m]*(2/3.);
	    }
	  //dF += -(dx*x*2*dens + drdR*x*x*ddens)*(xxdens[a1] + xxdens[a2])*2 - (dy*y*2*dens + drdR*y*y*ddens)*(yydens[a1] + yydens[a2])*2 - (dz*z*2*dens + drdR*z*z*ddens)*(zzdens[a1] + zzdens[a2])*2 - (dx*y*dens + dy*x*dens + drdR*x*y*ddens)*(xydens[a1] + xydens[a2])*4 - (dy*z*dens + dz*y*dens + drdR*y*z*ddens)*(yzdens[a1] + yzdens[a2])*4 - (dz*x*dens + dx*z*dens + drdR*z*x*ddens)*(zxdens[a1] + zxdens[a2])*4 + (drdR*ddens)*(dens2[a1] + dens2[a2])*(2/3.);
	}
      //for (int m = 0; m < n; m++)
      //	cerr << dFX[m][0] - dF[m][0] << " " << dFX[m][1] - dF[m][1] << " " << dFX[m][2] - dF[m][2] << endl; 
      for (int m = 0; m < n; m++)
	{
	  int a2 = neighbors[m];
	  for (int i = 0; i < 3; i++)
	    forces[a1][i] += dFX[m][i];
	  if (a2 < nAtoms)
	    for (int i = 0; i < 3; i++)
	      forces[a2][i] -= dFX[m][i];
	  if (counters.stresses != atoms->GetChangeCounter())
	    {
	      for (int i = 0; i < 3; i++)
		for (int j = i; j < 3; j++)
		  {
		    double ds = 0.5 * dFX[m][i] * diffs[m][j];
		    int k = stresscomp[i][j];
		    //stress[a1][k] += ds;
		    potstresses[a1][k] += ds;
		    if (a2 < nAtoms)
		      {
			//stress[a2][k] += ds;
			potstresses[a2][k] += ds;
		      }
		  }
	    }
	}
    }
  //  for (int a1 = 0; a1 < nAtoms; a1++)
  // for(int i = 0;i<6;i++)
  //   stress[a1][i] = potstresses[a1][i];	      
  counters.stresses = atoms->GetChangeCounter();
  counters.forces = atoms->GetChangeCounter();
  DEBUGPRINT;
}

void MoPotential::CalculateRhos()
  
{
  // cerr << "void MoPotential::CalculateRhos()" << endl;
  const int maxSize = 100;
  int neighbors[maxSize];
  Vec diffs[maxSize];
  double diffs2[maxSize];
  double r[maxSize];
  double invr[maxSize];
  double X[maxSize][3];
  Vec drdR[maxSize];
  double dE[maxSize];
  Vec dF[maxSize];
  double tmpdens0[maxSize];
  double funcvalue[maxSize];
  double tmpXdens[maxSize][3];
  double tmpXXdens[maxSize][6];
  memset(&atomicEnergies[0], 0, nAtoms * sizeof(double));
  memset(&potforces[0][0], 0, nAtoms * sizeof(Vec));
  double (*potstresses)[6];
  potstresses = (double (*)[6]) &(this->potstresses[0]);
  memset(&potstresses[0][0], 0, 6 * nAtoms * sizeof(double));
  memset(&lat[0], 0, nAtoms * sizeof(double));
  memset(&dens0[0], 0, nAtoms * sizeof(double));
  memset(&Factor[0], 0, nAtoms * sizeof(double));
  memset(&rho1[0], 0, nAtoms * sizeof(double));
  double (*Xdens)[3];
  Xdens  = (double (*)[3]) &(this->Xdens[0]);
  memset(&Xdens[0][0], 0, 3 * nAtoms * sizeof(double));
  memset(&rho2[0], 0, nAtoms * sizeof(double));
  memset(&dens2[0], 0, nAtoms * sizeof(double));
  double (*XXdens)[6];
  XXdens  = (double (*)[6]) &(this->XXdens[0]);
  memset(&XXdens[0][0], 0, 6 * nAtoms * sizeof(double));
  DEBUGPRINT;
  for (int a1 = 0; a1 < nAtoms; a1++)
    {
      int size = maxSize;
      int n = neighborList->
	GetNeighbors(a1, neighbors, diffs, diffs2, size);
      assert(size >= 0);
      for (int m = 0; m < n; m++)
	{
	  r[m] = sqrt(diffs2[m]);
	  invr[m] = 1/r[m];
	}
      for (int m = 0; m < n; m++)
	{
	  drdR[m]  = -diffs[m]*invr[m];
	  for (int i = 0; i < 3; i++)
	    X[m][i] = -drdR[m][i];
	}
      for (int m = 0; m < n; m++)
	{
	  double *pt_e, *pt_df;
	  //Vec drdR = -diffs[m] * invr[m];
	  double df;// = 2*dPairPotential(r);
	  pt_e = &dE[m];
	  pt_df = &df;
	  PairPotentialPointer(r[m], pt_e, pt_df);
	  //e = *pt_e;
	  //df = *pt_df*2;
	  dF[m] = -drdR[m] * df*2;	
	}
      for (int m = 0; m < n; m++)
	{
	  atomicEnergies[a1] += dE[m];
	  //forces[a1] += dF[m];
	  potforces[a1] += dF[m];
	  tmpdens0[m] = Density(r[m]);
	  funcvalue[m] = func(r[m]);
	}
      for (int m = 0; m < n; m++)
	for (int i = 0; i < 3; i++)
	  tmpXdens[m][i] = X[m][i]*funcvalue[m];	
      for (int m = 0; m < n; m++)
	for (int i = 0; i < 3; i++)
	  {
	    int j = (i+1)%3;
	    tmpXXdens[m][i] = X[m][i]*tmpXdens[m][i];
	    tmpXXdens[m][i+3] = X[m][j]*tmpXdens[m][i];
	  }
      for (int m = 0; m < n; m++)
	{
	  dens0[a1] += tmpdens0[m];
	  for (int i = 0; i < 3; i++)
	    Xdens[a1][i] += tmpXdens[m][i];
	  dens2[a1] += funcvalue[m];
	  for (int i = 0; i < 6; i++)
	    XXdens[a1][i] += tmpXXdens[m][i];
	}
      for (int m = 0; m < n; m++)
	{
	  int a2 = neighbors[m];
	  if (a2 < nAtoms)
	    {
	    atomicEnergies[a2] += dE[m];
	    //forces[a2] -= dF[m];
	    potforces[a2] -= dF[m];
	    dens0[a2] += tmpdens0[m];
	    for (int i = 0; i < 3; i++)
	      Xdens[a2][i] -= tmpXdens[m][i];
	    dens2[a2] += funcvalue[m];
	    for (int i = 0; i < 6; i++)
	      XXdens[a2][i] += tmpXXdens[m][i];
	    }
	  //if (stress)
	  //  {
	      for (int i = 0; i < 3; i++)
		for (int j = i; j < 3; j++)
		  {
		    double ds = 0.5 * dF[m][i] * diffs[m][j];
		    int k = stresscomp[i][j];
		    //stress[a1][k] += ds;
		    potstresses[a1][k] += ds;
		    if (a2 < nAtoms)
		      {
			//stress[a2][k] += ds;
			potstresses[a2][k] += ds;
		      }
		  }
	      //}
	}
    }
  DEBUGPRINT;
  for (int a1 = 0; a1 < nAtoms; a1++)
    {
      for (int i = 0; i < 3; i++)
	{
	  rho1[a1] += Xdens[a1][i]*Xdens[a1][i];// + ydens[a1]*ydens[a1] + zdens[a1]*zdens[a1];
	  rho2[a1] += XXdens[a1][i]*XXdens[a1][i] + 2*XXdens[a1][i+3]*XXdens[a1][i+3];
	}
      rho2[a1] -= dens2[a1]*dens2[a1]/3.;
      lat[a1] = LatticeConstant(dens0[a1]/norm);
//       if (a1 == 449360)
// 	cerr << dens0[a1]/norm<< " " << lat[a1] << " " << atomicEnergies[a1] << " " << EmbeddingEnergy(lat[a1]) << " " << PairPotentialRef(lat[a1]) << endl;
      atomicEnergies[a1] += EmbeddingEnergy(lat[a1]); 
      atomicEnergies[a1] -= PairPotentialRef(lat[a1]);
      Factor[a1] = (dEmbeddingEnergy(lat[a1])-dPairPotentialRef(lat[a1]))*dLatticeConstant(dens0[a1]/norm);
      atomicEnergies[a1] += c[0]*rho1[a1];
      atomicEnergies[a1] += c[1]*rho2[a1];
    }
  double maxlat=3.1;
  for (int a1 = 0; a1 < nAtoms; a1++)
    {
      if (lat[a1]>maxlat) 
	{
	  maxlat=lat[a1];
	}
    }  
  // cerr << "Maxlat" << maxlat << endl;
  //cerr << "densMin/latMax : "<< densMin << " / " << LatticeConstant(densMin) << " / " << Factor[iMin] << endl;
  //cerr << "densMax/latMin : "<< densMax << " / " << LatticeConstant(densMax) << " / " << Factor[iMax] << endl;
  DEBUGPRINT;
  GhostAtoms *ghostAtoms = dynamic_cast<GhostAtoms *>(atoms);
  if (ghostAtoms != 0)
    {
      // assert(dynamic_cast<ParallelAtoms *>(atoms) != 0);
  DEBUGPRINT;
  //  ghostAtoms->GetGhostPotential()->CommunicateData(&lat[0]);
  DEBUGPRINT;
      ghostAtoms->GetGhostPotential()->CommunicateData(&dens0[0]);
  DEBUGPRINT;
      ghostAtoms->GetGhostPotential()->CommunicateData(&Factor[0]);
  DEBUGPRINT;
  //ghostAtoms->GetGhostPotential()->CommunicateData(&rho1[0]);
  DEBUGPRINT;
      ghostAtoms->GetGhostPotential()->CommunicateData(&Xdens[0][0],3);
  DEBUGPRINT;
//       ghostAtoms->GetGhostPotential()->CommunicateData(&xdens[0]);
//       ghostAtoms->GetGhostPotential()->CommunicateData(&ydens[0]);
//       ghostAtoms->GetGhostPotential()->CommunicateData(&zdens[0]);
  //  ghostAtoms->GetGhostPotential()->CommunicateData(&rho2[0]);
  DEBUGPRINT;
      ghostAtoms->GetGhostPotential()->CommunicateData(&XXdens[0][0],6);
  DEBUGPRINT;
      ghostAtoms->GetGhostPotential()->CommunicateData(&dens2[0]);
//       ghostAtoms->GetGhostPotential()->CommunicateData(&xxdens[0]);
//       ghostAtoms->GetGhostPotential()->CommunicateData(&yydens[0]);
//       ghostAtoms->GetGhostPotential()->CommunicateData(&zzdens[0]);
//       ghostAtoms->GetGhostPotential()->CommunicateData(&xydens[0]);
//       ghostAtoms->GetGhostPotential()->CommunicateData(&yzdens[0]);
//       ghostAtoms->GetGhostPotential()->CommunicateData(&zxdens[0]);
      ghostAtoms->GetGhostPotential()->CommunicateData(&potforces[0][0],3);
      ghostAtoms->GetGhostPotential()->CommunicateData(&potstresses[0][0],6);
  DEBUGPRINT;
    }
  DEBUGPRINT;
}

