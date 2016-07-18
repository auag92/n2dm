/* bcuint.c - Copyright (c) 1998 Zyvex LLC.
* 
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions
* are met:
* 1. Redistributions of source code must retain the above copyright
*    notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in the
*    documentation and other materials provided with the distribution.
* 3. All advertising materials mentioning features or use of this software
*    or its derived works must display the following acknowledgement:
* 	This product includes software developed by Zyvex LLC.
* 
* This software is provided "as is" and any express or implied warranties,
* including, but not limited to, the implied warranties of merchantability
* or fitness for any particular purpose are disclaimed. In no event shall
* Zyvex LLC be liable for any direct, indirect, incidental, special,
* exemplary, or consequential damages (including, but not limited to,
* procurement of substitute goods or services; loss of use, data, or
* profits; or business interruption) however caused and on any theory of
* liability, whether in contract, strict liability, or tort (including
* negligence or otherwise) arising in any way out of the use of this
* software, even if advised of the possibility of such damage.
*/

#include "BrennerPotential.h"
#include "asapbrenner.h"
#include "Exception.h"
#include <math.h>
#include <stdio.h>
#include "expand.h"

/* Bicubic spline */

/* Hij(NiH,NiC)*/

#define MAX_BC_NEIGHBORS 32

static const char* ktype_name[5] = {NULL, "Carbon", "Hydrogen",
				    "Silicon", "Germanium"};

static int IN2[16+1][2+1];
/* Index is the inverse of IN2.
For II between 1 and 16 inclusive,
index[IN3[II][1],IN3[II][2]] == II. */
static int xindex [4][4];
static Double CLM[2+1][MAX_BC_NEIGHBORS+1][MAX_BC_NEIGHBORS+1][16+1];
static Float CLMS[2+1][MAX_BC_NEIGHBORS+1][MAX_BC_NEIGHBORS+1][16+1];

namespace ASAPSPACE {
extern int inter2d_in2_index[72][3];
extern double inter2d_in2_data[72][16];
}

void BrennerPotential::init_in2()
{
  DEBUGPRINT;
	int IC=0;
	int I, J, K, L, M;
	for (I=1; I<5; I++) {
		for (J=1; J<5; J++) {
			IC=IC+1;
			IN2[IC][1]=I-1;
			IN2[IC][2]=J-1;
			xindex[IN2[IC][1]][IN2[IC][2]]=IC;
		}
	}

	/* Zero bicubic spline coefficients*/
	for (I=1; I<3; I++)
		for (L=1; L <= MAX_BC_NEIGHBORS; L++)
			for (M=1; M <= MAX_BC_NEIGHBORS; M++)
				for (J=1; J<17; J++)
					CLM[I][L][M][J]=0.0;
	for (K=1; K<73; K++) {
                I = inter2d_in2_index[K-1][0];
                L = inter2d_in2_index[K-1][1];
                M = inter2d_in2_index[K-1][2];
		for (J=0; J<16; J++) {
			CLM[I][L][M][J+1] = (Float) inter2d_in2_data[K-1][J];
		}
	}
	{
		int littlei;
		for (littlei = 1; littlei < 3; littlei++)  {
			if(xindex[0][0]!=1) {
				int a = 10;
			}
			for (L=1; L <= MAX_BC_NEIGHBORS; L++)  {
				if(xindex[0][0]!=1) {
					int a = 10;
				}
				for (M=1; M <= MAX_BC_NEIGHBORS; M++) {
					int II;
					Double copy [16+1];
					if(xindex[0][0]!=1) {
						int a = 10;
					}
					for (II = 1; II < 17; II++) 
						copy [II] = 0.0;
					for (II = 1; II < 17; II++) {
						/* i and j are the exponents of XX1, XX2 for the
						term at position II in array CLM[L][M]. */
						int i, j;
						/* ti [z] is the coefficient of F1**z in the expansion of
						(L+F1)**i. Analogously for tj. */
						int ti[4], tj[4];
						/* ie will be the term of the expansion of (L+F1)**i that we'll
						be looking at.  Analogously je. */
						int ie, je, ke;
						i = IN2[II][1];
						j = IN2[II][2];
						expand (L, ti, i);
						expand (M, tj, j);
						if(xindex[0][0]!=1) {
		 					int a = 10;
						}
						for (ie = 0; ie < 4; ie++)
							for (je = 0; je < 4; je++) {
								copy [xindex [ie][je]] +=
									CLM[littlei][L][M][II]*ti[ie]*tj[je];
							}
					}
					if(xindex[0][0]!=1) {
						int a = 10;
			        }
                                        //DDL:
					//BUG? Why count to 65? although copy is only 17 bytes 
					//     AND data out of CLMS is accessed!!!
					for (II = 1; II < 17; II++) {
						CLMS [littlei][L][M][II] = copy [II];
						if(xindex[0][0]!=1) {
							int a = 10;
						}
					}
				}
			}
		}
	}
}

Float BrennerPotential::BCUINT(int KJ, Float XX1, Float XX2,
			       Float *ansy1_ptr, Float *ansy2_ptr)
{
	int NH, NC;
	/*bicubic spline*/

	if (XX1 < 1) XX1 = 1;
	if (XX2 < 1) XX2 = 1;
	NH=(int)floor(XX1);
	NC=(int)floor(XX2);

	assert (0 != NH);
	assert (0 != NC);
	if(NH > MAX_BC_NEIGHBORS || NC > MAX_BC_NEIGHBORS)
	{
	  throw AsapError("BrennerPotential: A ")
	    << (KJ >= 1 && KJ <= 4 ? ktype_name[KJ] : "invalid type")
	    << " atom has too many neighbors for the bicubic spline. It has "
	    << NH << " Hydrogen and " << NC
	    << " Carbon neighbors; the maximum is "
	    << MAX_BC_NEIGHBORS << "(MAX_BC_NEIGHBORS).";
	}
	if (KJ==0) {
	  throw AsapError("BrennerPotential: error BCUINT unexpected zero: KJ ")
	    << KJ << " NH " << NH << " NC " << NC
	    << " XX1 " << XX1 << " XX2 " << XX2;
	}
	{
		Float ANSY;
		Float ansy1 = 0;
		Float ansy2 = 0;
		const Float *clms_ptr = &CLMS[KJ][NH][NC][0];
		Float XX1F = XX1 - NH;
		Float XX2F = XX2 - NC;
		if (XX1F < 0.000001 && XX2F < 0.000001) {
			ANSY=clms_ptr [xindex[0][0]];
			ansy1=clms_ptr [xindex[1][0]];
			ansy2=clms_ptr [xindex[0][1]];
		} else {
			int j;
			Float xx1_pow[4];
			Float xx2_pow[4];
			Float xx1_d[4];
			Float xx2_d[4];
			xx1_pow[0] = xx2_pow[0] = 1;
			xx1_d[0] = xx2_d[0] = 0;
			for (j = 1; j < 4; j++) {
				xx1_pow[j] = xx1_pow[j-1]*XX1F;
				xx1_d[j] = xx1_pow[j-1]*j;
				xx2_pow[j] = xx2_pow[j-1]*XX2F;
				xx2_d[j] = xx2_pow[j-1]*j;
			}
			ANSY=0;
			for (j=1; j<17; j++) {
				const Float c = *++clms_ptr;
				ANSY += c*xx1_pow[IN2[j][1]]*xx2_pow[IN2[j][2]];
				ansy1 += c*xx1_d[IN2[j][1]]*xx2_pow[IN2[j][2]];
				ansy2 += c*xx1_pow[IN2[j][1]]*xx2_d[IN2[j][2]];
			}
		}
		if(ansy1_ptr) *ansy1_ptr = ansy1;
		if(ansy2_ptr) *ansy2_ptr = ansy2;
		return ANSY;
	}
}

/* Compute P sub i j, which is referenced in equation 7 on page 10.  
KJ is probably the identity of the "j" atom in P sub i j.  I don't know why
we don't have a KI argument.
Judging by the variable names, XNT1 is probably N super H sub i and XNT2 is
N super C sub i.  Note that this differs from the order of arguments of P
sub i j in the paper.  
ansy1_ptr is an output, the derivative of P sub i j with respect to
XNT1.  
ansy2_ptr is an output, the derivative of P sub i j with respect to
XNT2.  */

