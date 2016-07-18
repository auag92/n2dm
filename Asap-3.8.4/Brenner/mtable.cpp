/*
 This code was released into the public domain by Peter McCluskey on 6/2/2000.
 It is adapted from code released into the public domain by Donald Brenner on
 2/29/2000.
*/

#include "BrennerPotential.h"
#include "asapbrenner.h"
#include <math.h>

/*
 * All array subscripts in this routine are 1 less than in the equivalent
 * Fortran code
*/

void BrennerPotential::mtable(Float drtable[4][4][NTAB],
			      Float ddtab[4][4],
			      Float rtable[4][4][NTAB],
			      Float atable[4][4][NTAB],
			      Float datable[4][4][NTAB], 
			      Float tabfc[4][4][NTAB],
			      Float tabdfc[4][4][NTAB])
{
  /* generate lookup tables for bond-order potentials */
  /* This is all only done once, so use double precision for all temporary
     variables. */
  int ki;
  int i, j;
  Double va, vb, vc;
  Double dva, dvb, dvc;
  Double vv, dvv, dvm;
  Double ff1, ff2, df1, df2;
  
  /* PID [i][j] will be pi / ((D super max sub i j) - (D super min sub i j)),
     which is used in computing f super c sub i j (r) in equation 18 on page
     16. */
  Double PID[4][4];
  /* AD is B sub i j 1, as in equation 5 on page 7. */
  Double AD[4][4];
  /* AXL is Beta sub i j 1, as in equation 5 on page 7. */
  Double AXL[4][4];
  /* BD is B sub i j 2, as in equation 5 on page 7. */
  Double BD[4][4];
  /* BXL is beta sub i j 2, as in equation 5 on page 7. */
  Double BXL[4][4];
  /* CD is B sub i j 2, as in equation 5 on page 7. */
  Double CD[4][4];
  /* CXL is beta sub i j 2, as in equation 5 on page 7. */
  Double CXL[4][4];
  Double DD[4][4];
  Double DXL[4][4];
  Double ED[4][4];
  /* RB1[i][j] is D super min sub i j in equation 18 on page 16.  When we're
     computing f super c sub ij, and the radius is less than 
     D super min sub i j, then the f super c sub i j is 1. */
  Double RB1[4][4];
  Double CHI[4][4];
  
  for(i = 0; i < 4; ++i)
    for(j = 0; j < 4; ++j)
      PID[i][j] = AD[i][j] = AXL[i][j] = BD[i][j] = BXL[i][j]
        = CD[i][j] = CXL[i][j] = DD[i][j] = DXL[i][j] = ED[i][j] = RB1[i][j]
        = CHI[i][j] = 0.0;
  
  /* Carbon */
  
  AD[0][0] = 12388.79197798375;
  AXL[0][0] = 4.720452312717397;
  BD[0][0] = 17.56740646508968;
  BXL[0][0] = 1.433213249951261;
  CD[0][0] = 30.71493208065162;
  CXL[0][0] = 1.382691250599169;
  DD[0][0] = 10953.54416216992;
  DXL[0][0] = 4.746539060659529;
  ED[0][0] = 0.3134602960832605;
  RB1[0][0] = 1.7;
  rb2[0][0] = 2.0;
  PID[0][0] = M_PI/(rb2[0][0] -RB1[0][0]);
  
  /* Hydrogen */
  
  AD[1][1] = 29.6325931;
  AXL[1][1] = 1.715892169856421;
  BD[1][1] = 0.0;
  BXL[1][1] = 1.0;
  CD[1][1] = 0.0;
  CXL[1][1] = 1.0;
  DD[1][1] = 32.81735574722296;
  DXL[1][1] = 3.536298648376465;
  ED[1][1] = 0.3704714870452888;
  RB1[1][1] = 1.10;
  rb2[1][1] = 1.70;  
  PID[1][1] = M_PI/(rb2[1][1]-RB1[1][1]);
  
  /* CARBON-HYDROGEN */
  
  AD[1][0] = 32.35518665873256;
  AXL[1][0] = 1.434458059249837;
  DD[1][0] = 149.9409872288120;
  DXL[1][0] =  4.102549828548784;
  ED[1][0] = 0.3407757282257080;
  
  BD[1][0] = 0.0;
  BXL[1][0] = 1.0;
  CD[1][0] = 0.0;
  CXL[1][0] = 1.0;
  
  AD[0][1] = AD[1][0];
  AXL[0][1] = AXL[1][0];
  BD[0][1] = BD[1][0];
  BXL[0][1] = BXL[1][0];
  CD[0][1] = CD[1][0];
  CXL[0][1] = CXL[1][0];
  DD[0][1] = DD[1][0];
  DXL[0][1] = DXL[1][0];
  ED[0][1] = ED[1][0];
  
  RB1[1][0] = 1.3;
  rb2[1][0] = 1.8;
  PID[1][0] = M_PI/(rb2[1][0] - RB1[1][0]);
  /* Note that PIDT really is PI/0.3 in pibond.c. */
  /*	PIDT=PI/0.3; */
  
  RB1[0][1] = RB1[1][0];
  rb2[0][1] = rb2[1][0];
  PID[0][1] = M_PI/(rb2[0][1] - RB1[0][1]);
  /*
   *  TERSOFF-III SILICON
   */
  DXL[2][2] = 2.4799;
  AXL[2][2] = 1.7322;
  DD[2][2] = 1830.8;
  AD[2][2] = 471.18;
  
  RB1[2][2] = 2.7;
  rb2[2][2] = 3.0;
  PID[2][2] = M_PI/(rb2[2][2] - RB1[2][2]);
  /*
   *  TERSOFF GERMANIUM
   */
  DXL[3][3] = 2.4451;
  AXL[3][3] = 1.7047;
  DD[3][3] = 1769.0;
  AD[3][3] = 0.50 * 419.23;
  
  CHI[2][2] = 1.0;
  PID[3][3] = M_PI/(rb2[3][3] - RB1[3][3]);
  CHI[3][3] = 1.0;
  CHI[3][2] = 1.00061;
  CHI[2][3] = CHI[3][2];
  
  RB1[3][3] = 2.7;
  rb2[3][3] = 3.0;
  
  /*
   * Mixed SILICON-GERMANIUM
   */
  DXL[3][2] = (DXL[3][3]+DXL[2][2])/2.0;
  DXL[2][3] = DXL[3][2];
  AXL[3][2] = (AXL[3][3]+AXL[2][2])/2.0;
  AXL[2][3] = AXL[3][2];
  DD[3][2] = sqrt(DD[3][3]*DD[2][2]);
  DD[2][3] = DD[3][2];
  AD[3][2] = sqrt(AD[2][2]*AD[3][3]);
  AD[2][3] = AD[3][2];
  
  RB1[3][2] = sqrt(RB1[2][2]*RB1[3][3]);
  RB1[2][3] = RB1[3][2];
  rb2[3][2] = sqrt(rb2[2][2]*rb2[3][3]);
  rb2[2][3] = rb2[3][2];
  PID[3][2] = M_PI/(rb2[3][2]-RB1[3][2]);
  PID[2][3] = PID[3][2];
  
  
  for(ki = 0; ki < 4; ++ki) {
    int kj;
    for(kj = ki; kj < 4; ++kj) {
      for(i = 1; i < NTAB-1; ++i) {
        if(rb2[ki][kj] != 0) {
          /* cut-off function.  That is, f super c sub i j (r) defined in
             equation 18 on page 16. */
          Double FC=0.0;
          /* Derivative of the cut-off function with respect to r. */
          Double DFC=0.0;
          Double rc, rsq;
          ddtab[ki][kj] = (NTAB - 2)/rb2[ki][kj];
          ddtab[kj][ki] = ddtab[ki][kj];
          /* This used to compute rc by adding ddtab[ki][kj] to the previous
             value of rc, but I figured that was bad because it lets roundoff
             errors accumulate, so now we multiply instead.  Tim Freeman 5
             Aug 2000. */
          rc = i / ddtab[ki][kj];
          rsq = rc*rc;
          
          if(rc < rb2[ki][kj]) {
            Double DTEMP = PID[ki][kj]*(rc-RB1[ki][kj]);
            FC = (1.0 + cos(DTEMP))/2.0;
            DFC = -PID[ki][kj]/2.0*sin(DTEMP);
          }
          
          if(rc <= RB1[ki][kj]) {
            FC=1.0;
            DFC=0.0;
          }
          
          tabfc[ki][kj][i] = FC;
          tabfc[kj][ki][i] = tabfc[ki][kj][i];
          tabdfc[ki][kj][i] = DFC;
          tabdfc[kj][ki][i] = tabdfc[ki][kj][i];
          /* attractive pair terms */
          va = AD[ki][kj]*exp(-AXL[ki][kj]*rc);
          /* dva is partial of va with respect to rc. */
          dva = -AXL[ki][kj]*va;
          
          vb = BD[ki][kj]*exp(-BXL[ki][kj]*rc);
          dvb =- BXL[ki][kj]*vb;
          
          vc = CD[ki][kj]*exp(-CXL[ki][kj]*rc);
          dvc = -CXL[ki][kj]*vc;
          
          /* vv is the sigma in equation 5 page 7, except I don't know why
             we're dividing by 2. */
          vv = (va+vb+vc)/2.0;
          /* dvv is the partial of vv with respect to rc. */
          dvv = (dva+dvb+dvc)/2.0;
          /* atable has values for V super A sub i j. */
          atable[ki][kj][i] = FC*vv;
          atable[kj][ki][i] = atable[ki][kj][i];
          /* datable has values for the partial V super A sub i j wrt rc. */
          datable[ki][kj][i] = (FC*dvv+DFC*vv)/rc;
          datable[kj][ki][i] = datable[ki][kj][i];
          /* repulsive pair terms */
          ff1 = DD[ki][kj]*exp(-DXL[ki][kj]*rc);
          df1 = -DXL[ki][kj]*ff1;
          
          ff2 = (1.0+ED[ki][kj]/rc);
          df2 = -ED[ki][kj]/rsq;
          
          vv = ff1*ff2;
          dvm = (df1*ff2 + ff1*df2);
          /* Next line is equation 4 on page 7. */
          rtable[ki][kj][i] = vv*FC;
          rtable[kj][ki][i] = rtable[ki][kj][i];
          drtable[ki][kj][i] = -(FC*dvm+DFC*vv)/rc;
          drtable[kj][ki][i] = drtable[ki][kj][i];
        } else { /* rb2 is 0. */
          tabfc[ki][kj][i] = 0.0;
          tabfc[kj][ki][i] = tabfc[ki][kj][i];
          tabdfc[ki][kj][i] = 0.0;
          tabdfc[kj][ki][i] = tabdfc[ki][kj][i];
          atable[ki][kj][i] = 0.0;
          atable[kj][ki][i] = atable[ki][kj][i];
          datable[ki][kj][i] = 0.0;
          datable[kj][ki][i] = datable[ki][kj][i];
          rtable[ki][kj][i] = 0.0;
          rtable[kj][ki][i] = rtable[ki][kj][i];
          drtable[ki][kj][i] = 0.0;
          drtable[kj][ki][i] = drtable[ki][kj][i];
        }
      } /* End of the for loop over i. */
      atable[ki][kj][0] = atable[ki][kj][1];
      atable[kj][ki][0] = atable[ki][kj][0];
      datable[ki][kj][0] = datable[ki][kj][1];
      datable[kj][ki][0] = datable[ki][kj][0];
      rtable[ki][kj][0] = rtable[ki][kj][1];
      rtable[kj][ki][0] = rtable[ki][kj][0];
      drtable[ki][kj][0] = drtable[ki][kj][1];
      drtable[kj][ki][0] = drtable[ki][kj][0];
      tabfc[ki][kj][0] = 0.0;
      tabfc[kj][ki][0] = 0.0;
      tabdfc[ki][kj][0] = 0.0;
      tabdfc[kj][ki][0] = 0.0;
      
      atable[ki][kj][NTAB-1] = 0.0;
      atable[kj][ki][NTAB-1] = atable[ki][kj][NTAB-1];
      datable[ki][kj][NTAB-1] = 0.0;
      datable[kj][ki][NTAB-1] = datable[ki][kj][NTAB-1];
      rtable[ki][kj][NTAB-1] = 0.0;
      rtable[kj][ki][NTAB-1] = rtable[ki][kj][NTAB-1];
      drtable[ki][kj][NTAB-1] = 0.0;
      drtable[kj][ki][NTAB-1] = drtable[ki][kj][NTAB-1];
      tabfc[ki][kj][NTAB-1] = 0.0;
      tabfc[kj][ki][NTAB-1] = tabfc[ki][kj][NTAB-1];
      tabdfc[ki][kj][NTAB-1] = 0.0;
      tabdfc[kj][ki][NTAB-1] = tabdfc[ki][kj][NTAB-1];
    } /* End loop over kj. */
  } /* End loop over ki. */
} /* End definition of mtable. */
