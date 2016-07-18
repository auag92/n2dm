// -*- C++ -*-
// Make Radial Distribution Functions (RDFs) grouped by elements and
// possibly by another integer tag, for example encoding the location.
// Only the counting is done here.  The normalization etc is done in Python.
//
// Copyright (C) 2007-2011 Jakob Schiotz and Center for Individual
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

#ifndef _RAWRADIALDISTRIBUTION_H
#define _RAWRADIALDISTRIBUTION_H

#include "Asap.h"
#include <vector>
#include <map>
#include <utility>
using std::vector;
using std::map;
using std::pair;

namespace ASAPSPACE {

class Atoms;

typedef vector< map< pair<int,int>, long* > > RDFtype;
typedef vector< map< int, long > > RDFcountType;

/// Produce an RDF but leave the normalization to a Python object.

/// This function does the counting needed to produce radial
/// distribution functions.  Atoms may be divided into groups, then
/// separate RDFS are produced for each group.  Partial RDFs for each
/// pair of elements are also calculated.
///
/// Input parameters:
///   atoms:  The atoms object to be analyzed.
///   ngroups: The number of groups.
///   groups: An interger for each atom, giving the group.  If no groups
///           are wanted, place all atoms in group 0.
///   nBins: Number of bins in the histogram.
///
/// Output parameters:
///   result: Contains the RDF for atom pairs (A,B) in group G in
///           result[G][pair(A,B)]
///   counts: Contains the number of atoms of element A in group G in
///           counts[G][A]
///   globalresult: The global RDF.
void RawRadialDistribution(PyObject *pyatoms, int nGroups, const int *group,
			   double rCut, int nBins,
			   RDFtype &result, RDFcountType &counts,
			   long *globalresult);

} // end namespace

#endif // _RAWRADIALDISTRIBUTION_H
