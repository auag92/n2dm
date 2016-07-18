// Copyright (C) 2008-2011 Jakob Schiotz and Center for Individual
// Nanoparticle Functionality, Department of Physics, Technical
// University of Denmark.  Email: schiotz@fysik.dtu.dk
//
// This file is part of Asap version 3.
// Asap is released under the GNU Lesser Public License (LGPL) version 3.
// However, the parts of Asap distributed within the OpenKIM project
// (including this file) are also released under the Common Development
// and Distribution License (CDDL) version 1.0.
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


#include "NeighborList.h"
#include "Atoms.h"
#include "Exception.h"
#include "Timing.h"
//#define ASAPDEBUG
#include "Debug.h"
#include <iostream>
#include <stdio.h>
using std::cerr;
using std::endl;
using std::flush;

extern int verbose;

//#define MAXLIST 500
//#define CHKFULLLIST
//#define CHECKCONSISTENCY
//#define PRINTLISTSIZE

NeighborList::NeighborList(Atoms *a, double rCut, double driftfactor)
{
  CONSTRUCTOR;
  if (a == NULL)
    atoms = new NormalAtoms();
  else
    {
      atoms = a;
      AsapAtoms_INCREF(atoms);
    }
  nAtoms = 0;  // Not valid yet
  invalid = true;
  firsttime = true;
  fulllists = false;
  this->rCut = rCut;
  rCut2 = rCut * rCut;
  drift = driftfactor * rCut;
  drift2 = drift * drift;
  // A cell list is used to build the neighbor list.  Its cutoff must
  // include the drift of this list, but it does not need its own
  // drift.  It is marked as a slave, so it does not trigger
  // additional migrations when it updates.
  if (rCut <= 0.0)
    throw AsapError("NeighborList: cutoff distance must be positive.");
  PyAsap_NeighborLocatorObject *nbl;
  nbl = PyAsap_NewNeighborCellLocator(atoms, rCut+2*drift, 0.0, true);
  cells_obj = (PyObject *) nbl;
  cells = dynamic_cast<NeighborCellLocator*>(nbl->cobj);
  assert(cells);
  cells->GetTranslationTable(translationTable);
}

NeighborList::~NeighborList()
{
  DESTRUCTOR;
  CHECKREF(cells_obj);
  Py_DECREF(cells_obj);
  AsapAtoms_DECREF(atoms);
}

void NeighborList::EnableFullNeighborLists()
{
  invalid = true;
  fulllists = true;
}

void NeighborList::MakeList()
{
  RETURNIFASAPERROR;
  USETIMER("NeighborList::MakeList");
  // Do we have any ghosts?
  //GhostAtoms *ghostAtoms = dynamic_cast<GhostAtoms *>(atoms);

  const Vec *ss = NULL;
  const bool *periodic = NULL;
  const double *superCellHeights = NULL;
  DEBUGPRINT;

#ifdef _OPENMP
#pragma omp single copyprivate(ss,periodic,superCellHeights)
#endif // _OPENMP
  {
    ss = atoms->GetCell();
    periodic = atoms->GetBoundaryConditions();
    superCellHeights = atoms->GetCellHeights();

    DEBUGPRINT;
    if (verbose >= 1)
      cerr << " NeighborList-Update ";
    nAtoms = atoms->GetNumberOfAtoms();
    nAllAtoms = nAtoms + atoms->GetNumberOfGhostAtoms();
    if (nAllAtoms != nbList.size())
      nbList.resize(nAllAtoms);
    if (fulllists && (nAllAtoms != complNbList.size()))
      complNbList.resize(nAllAtoms);
    maxLength = 0;
    memcpy(storedSuperCell, ss, 3*sizeof(Vec));
    DEBUGPRINT;
    for (int i = 0; i < 3; i++)
      {
        pbc[i] = periodic[i];     // Store for later tests
        referenceSuperCell[i] = ss[i];
      }
    // Check the size.
    for (int i = 0; i < 3; i++)
      if (periodic[i] && superCellHeights[i] < 2 * rCut)
        THROW( AsapError("The height of the cell (")
        << superCellHeights[i] << ") must be larger than " << 2 * rCut );

    DEBUGPRINT;
    cells->UpdateNeighborList();
    DEBUGPRINT
    assert(nAtoms == atoms->GetNumberOfAtoms());
  }
  RETURNIFASAPERROR;
  
  // Make the list
  int myMaxLength = 0;
#ifdef PRINTLISTSIZE
  int totlistsize = 0;
#endif
  if (firsttime)
    {
      vector< pair<int, translationsidx_t> > buf;
      buf.reserve(cells->MaxNeighborListLength());
#ifdef _OPENMP
#pragma omp for nowait
#endif // _OPENMP
      for (int i = 0; i < nAllAtoms; i++)
	{
	  int l = cells->GetListAndTranslations(i, buf);
#ifdef MAXLIST
	  if (l > MAXLIST)
	    THROW( AsapError("Unreasonably long neighbor list for atom ")
	      << i << ": " << l << " elements." );
#endif
	  nbList[i].reserve(l + l/20 + 2);
	  nbList[i].clear();
	  nbList[i].insert(nbList[i].begin(), buf.begin(), buf.end());
	  if (fulllists)
	    {
	      int l2 = cells->GetComplementaryListAndTranslations(i, buf);
	      l += l2;
	      complNbList[i].reserve(l2 + l2/20 + 2);
	      complNbList[i].clear();
	      complNbList[i].insert(complNbList[i].begin(),
				    buf.begin(), buf.end());
	    }
	  if (l > myMaxLength)
	    myMaxLength = l;
#ifdef PRINTLISTSIZE
	  totlistsize += l;
#endif
	}
    }
  else
    {
#ifdef _OPENMP
#pragma omp for nowait
#endif // _OPENMP
      for (int i = 0; i < nAllAtoms; i++)
	{
	  int l = cells->GetListAndTranslations(i, nbList[i]);
#ifdef MAXLIST
	  if (l > MAXLIST)
	    THROW( AsapError("Unreasonably long neighbor list for atom ")
	      << i << ": " << l << " elements." );
#endif
	  if (fulllists)
	    l += cells->GetComplementaryListAndTranslations(i, complNbList[i]);
	  if (l > myMaxLength)
	    myMaxLength = l;
#ifdef PRINTLISTSIZE
          totlistsize += l;
#endif
	}
    }
  myMaxLength++;   // We need space in arrays for a candidate atom after the last neighbor.
#ifdef _OPENMP
#pragma omp critical
#endif // _OPENMP
  {
    if (myMaxLength > maxLength)
      maxLength = myMaxLength;
  }
#ifdef _OPENMP
#pragma omp single
#endif // _OPENMP
  {
    DEBUGPRINT;
    invalid = false;
    firsttime = false;
  }
#ifdef CHKFULLLIST
  CheckFullListConsistency("MakeList");
#endif
#ifdef PRINTLISTSIZE
  std::cerr << "NeighborList::MakeList: nAtoms=" << nAtoms << "  nAllAtoms=" << nAllAtoms << " MaxLength=" << maxLength <<  "  ListTotal=" << totlistsize << std::endl;
#endif
  DEBUGPRINT;
}

bool NeighborList::CheckNeighborList()
{
  RETURNIFASAPERROR2(false);
  USETIMER("NeighborList::CheckNeighborList");

  if (invalid)
    return true;

  int n_at = atoms->GetNumberOfAtoms();
  if (nAtoms != n_at || nAllAtoms != n_at + atoms->GetNumberOfGhostAtoms())
    return true;

  bool updateRequired = false;
  DEBUGPRINT;
  cells->RenormalizePositions(); // XXX Move out of omp single once parallelized.
  DEBUGPRINT;

  const bool *newpbc = atoms->GetBoundaryConditions();

  if (invalid && verbose)
    cerr << "NeighborList::CheckAndUpdateNeighborList: NBList has been marked invalid." << endl;

  const Vec *ss = atoms->GetCell();
  memcpy(storedSuperCell, ss, 3*sizeof(Vec));

  if (nAtoms != atoms->GetNumberOfAtoms() || pbc[0] != newpbc[0]
     || pbc[1] != newpbc[1] || pbc[2] != newpbc[2])
    {
      invalid = true;
    }

  if (invalid)
    cells->Invalidate();

  updateRequired = invalid;
  const Vec *positions = atoms->GetPositions();
  const Vec *referencePositions = cells->GetReferencePositions();

  if (!updateRequired)
    {
      // Check how much the cell has changed.
      double max_strain_disp = GetMaxStrainDisplacement();
      double max2 = drift - max_strain_disp;
      if (max2 <= 0.0)
        updateRequired = true;
      else
        {
          max2 = max2 * max2;
          for (int n = 0; n < nAtoms; n++)
            if (Length2(positions[n] - referencePositions[n]) > max2)
              {
                updateRequired = true;
                break;
              }
        }
    }

#ifdef CHKFULLLIST
if (!updateRequired)
  CheckFullListConsistency("CheckNeighborList");
#endif
return updateRequired;
}

void NeighborList::UpdateNeighborList()
{
  MakeList();
}

bool NeighborList::CheckAndUpdateNeighborList()
{
  RETURNIFASAPERROR2(false);
  DEBUGPRINT;
  bool update = CheckNeighborList();
  MEMORY;
  if (update)
    UpdateNeighborList();
  return update;
}

bool NeighborList::CheckAndUpdateNeighborList(PyObject *atoms_obj)
{
  atoms->Begin(atoms_obj);
  CHECKNOASAPERROR;
  bool res = CheckAndUpdateNeighborList();
  PROPAGATEASAPERROR;
  atoms->End();
  return res;
}

// Some code could be saved by combining GetNeighbors and
// GetFullNeighbors but it would mean a slight performance cost in
// GetNeighbors, which is critical for overall performance.
int NeighborList::GetNeighbors(int a1, int *neighbors, Vec *diffs,
				double *diffs2, int& size, double r) const
{
  if (invalid)
    {
      DEBUGPRINT;
      THROW( AsapError("NeighborList has been invalidated, possibly by another NeighborList using the same atoms.") );
    }

  if (size < nbList[a1].size())
    {
      DEBUGPRINT;
      THROW( AsapError("NeighborList::GetNeighbors: list overflow.") );
    }
  RETURNIFASAPERROR2(0);

  const vector<Vec> &positions = cells->GetWrappedPositions();
  // Need to use GET_CELL instead of GetCell as the atoms are not open
  // when called from the Python interface.
  const Vec *superCell = atoms->GET_CELL();
  double rC2 = rCut2;
  if (r > 0.0)
    rC2 = r * r;
  Vec pos1 = positions[a1];
  int nNeighbors = 0;

  typedef vector< pair<int,translationsidx_t> >::const_iterator iterat;
  iterat terminate = nbList[a1].end();
  for (iterat a2 = nbList[a1].begin(); a2 < terminate; ++a2)
    {
      assert(a2->second >= 0 && a2->second < translationTable.size());
      assert(a2->first >= 0 && a2->first < nAllAtoms);
      assert(nNeighbors < size);
      const IVec &celltranslation = translationTable[a2->second];
      diffs[nNeighbors] = positions[a2->first] - pos1
	- celltranslation[0] * superCell[0]
	- celltranslation[1] * superCell[1]
	- celltranslation[2] * superCell[2];
      double d2 = Length2(diffs[nNeighbors]);
      if (d2 < rC2)
	{
	  diffs2[nNeighbors] = d2;
	  neighbors[nNeighbors] = a2->first;
	  nNeighbors++;
	}
    }
  size -= nNeighbors;
  assert(size >= 0);
  return nNeighbors;
}

void NeighborList::GetNeighbors(int a1, vector<int> &neighbors) const
{
  if (invalid)
    THROW( AsapError("NeighborList has been invalidated, possibly by another NeighborList using the same atoms.") );
  RETURNIFASAPERROR;

  neighbors.clear();
  const vector<Vec> &positions = cells->GetWrappedPositions();
  const Vec *superCell = atoms->GET_CELL();  // atoms may not be open
  Vec pos1 = positions[a1];

  typedef vector< pair<int,translationsidx_t> >::const_iterator iterat;
  iterat terminate = nbList[a1].end();
  for (iterat a2 = nbList[a1].begin(); a2 < terminate; ++a2)
    {
      const IVec &celltranslation = translationTable[a2->second];
      Vec diff = positions[a2->first] - pos1
	- celltranslation[0] * superCell[0]
	- celltranslation[1] * superCell[1]
	- celltranslation[2] * superCell[2];
      double d2 = Length2(diff);
      if (d2 < rCut2)
	  neighbors.push_back(a2->first);
    }
}

int NeighborList::GetFullNeighbors(int a1, int *neighbors, Vec *diffs,
				    double *diffs2, int& size, double r) const
{
  if (!fulllists)
    THROW( AsapError("Calling NeighborList::GetFullNeighbors but full lists are not enabled.") );
  
  if (invalid)
    THROW( AsapError("NeighborList has been invalidated, possibly by another NeighborList using the same atoms.") );

  if (size < nbList[a1].size() + complNbList[a1].size())
    {
      std::cerr << "Neighborlist OVERFLOW: a=" << a1 << " size=" << size << " len(nbl)=" << nbList[a1].size()
          << " len(cmpnbl)=" << complNbList[a1].size() << std::endl;
      THROW( AsapError("NeighborList::GetFullNeighbors: list overflow.") );
    }
  RETURNIFASAPERROR2(0);

  const vector<Vec> &positions = cells->GetWrappedPositions();
  // Need to use GET_CELL instead of GetCell as the atoms are not open
  // when called from the Python interface.
  const Vec *superCell = atoms->GET_CELL();
  double rC2 = rCut2;
  if (r > 0.0)
    rC2 = r * r;
  Vec pos1 = positions[a1];
  int nNeighbors = 0;

  typedef vector< pair<int,translationsidx_t> >::const_iterator iterat;
  iterat terminate = nbList[a1].end();
  for (iterat a2 = nbList[a1].begin(); a2 < terminate; ++a2)
    {
      const IVec &celltranslation = translationTable[a2->second];
      diffs[nNeighbors] = positions[a2->first] - pos1
	- celltranslation[0] * superCell[0]
	- celltranslation[1] * superCell[1]
	- celltranslation[2] * superCell[2];
      double d2 = Length2(diffs[nNeighbors]);
      if (d2 < rC2)
	{
	  diffs2[nNeighbors] = d2;
	  neighbors[nNeighbors] = a2->first;
	  nNeighbors++;
	}
    }
  terminate = complNbList[a1].end();
  for (iterat a2 = complNbList[a1].begin(); a2 < terminate; ++a2)
    {
      const IVec &celltranslation = translationTable[a2->second];
      diffs[nNeighbors] = positions[a2->first] - pos1
	- celltranslation[0] * superCell[0]
	- celltranslation[1] * superCell[1]
	- celltranslation[2] * superCell[2];
      double d2 = Length2(diffs[nNeighbors]);
      if (d2 < rC2)
	{
	  diffs2[nNeighbors] = d2;
	  neighbors[nNeighbors] = a2->first;
	  nNeighbors++;
	}
    }
  size -= nNeighbors;
  assert(size >= 0);
  return nNeighbors;
}

void NeighborList::GetFullNeighbors(int a1, vector<int> &neighbors) const
{
  if (!fulllists)
    THROW( AsapError("Calling NeighborList::GetFullNeighbors but full lists are not enabled.") );
  RETURNIFASAPERROR;
  
  const vector<Vec> &positions = cells->GetWrappedPositions();
  const Vec *superCell;
  if (atoms->IsActive())
    superCell = atoms->GetCell();
  else
    superCell = storedSuperCell;
  Vec pos1 = positions[a1];

  // Get the first half of the neighbor list
  GetNeighbors(a1, neighbors);

  // Get the second half
  typedef vector< pair<int,translationsidx_t> >::const_iterator iterat;
  iterat terminate = complNbList[a1].end();
  for (iterat a2 = complNbList[a1].begin(); a2 < terminate; ++a2)
    {
      const IVec &celltranslation = translationTable[a2->second];
      Vec diff = positions[a2->first] - pos1
	- celltranslation[0] * superCell[0]
	- celltranslation[1] * superCell[1]
	- celltranslation[2] * superCell[2];
      double d2 = Length2(diff);
      if (d2 < rCut2)
	  neighbors.push_back(a2->first);
    }
}

void NeighborList::RemakeLists(const set<int> &modified, set<int> &affected)
{
  // Never called in OpenMP context - throwing exceptions is OK.
#ifdef CHKFULLLIST
  CheckFullListConsistency("Beginning of RemakeLists", false);
#endif // CHKFULLLIST
#ifdef CHECKCONSISTENCY
  const vector<Vec> &pos = cells->GetWrappedPositions();
  const Vec *superCell = atoms->GetCell();
#endif // CHECKCONSISTENCY
  
  typedef vector< pair<int,translationsidx_t> >::iterator nbiterat;

  if (invalid)
    throw AsapError("NeighborList has been invalidated, possibly by another neighbor list sharing the same atoms.");
  if (!fulllists)
    throw AsapError("NeighborList::RemakeLists() only supported if 'full neigbor lists' are enabled.");
  
  assert(modified.size() > 0);
  affected.clear();
  affected.insert(modified.begin(), modified.end());

  // Update the NeighborCellLocator.  The positions of the
  // modified atoms are normalized here.
  cells->RemakeLists_Simple(modified);

  // Collect the current neighbors to the modified atoms, and
  // remove the atoms from the neighbors' lists.  Atoms that are
  // already in the modified set are not affected.
  for (set<int>::const_iterator i = modified.begin(); i != modified.end(); ++i)
    {
      nbiterat terminate = nbList[*i].end();
      for (nbiterat j = nbList[*i].begin(); j != terminate; ++j)
	{
	  // j->first is a neighbor to modified atom *i
	  affected.insert(j->first);
	  if (!modified.count(j->first))
	    {
	      // Not in modified set.  Fix its complementary list:
	      // Remove atom *i from compl. neighbor list of j->first 
	      //int x = complNbList[j->first].erase(*i);
	      nbiterat k = complNbList[j->first].begin();
	      while ((k->first != *i) && (k != complNbList[j->first].end()))
		++k;
	      if (k != complNbList[j->first].end())
		  complNbList[j->first].erase(k);
#ifdef CHECKCONSISTENCY
	      else
		{
		  // Looks like an inconsistency, but we have to live
		  // with them if the distance is above rCut.
		  IVec celltranslation = translationTable[j->second];
		  Vec pos1 = pos[*i] + celltranslation[0] * superCell[0]
		    + celltranslation[1] * superCell[1]
		    + celltranslation[2] * superCell[2];
		  if (Length2(pos1 - pos[j->first]) < rCut2)
		    throw AsapError("Inconsistent neighbor list data: Did not find ")
		      << *i << " on complementary nb list of atom "
		      << j->first;
		}
#endif // CHECKCONSISTENCY
	    }
	}
      terminate = complNbList[*i].end();
      for (nbiterat j = complNbList[*i].begin(); j != terminate; ++j)
	{
	  affected.insert(j->first);
	  if (!modified.count(j->first))
	    {
	      // Not in modified set.  Fix its complementary list
	      //int x = nbList[j->first].erase(*i);
	      nbiterat k = nbList[j->first].begin();
	      while ((k->first != *i) && (k != nbList[j->first].end()))
		++k;
	      if (k != nbList[j->first].end())
		  nbList[j->first].erase(k);
#ifdef CHECKCONSISTENCY
	      else
		{
		  // Looks like an inconsistency, but we have to live
		  // with them if the distance is above rCut.
		  IVec celltranslation = translationTable[j->second];
		  Vec pos1 = pos[*i] + celltranslation[0] * superCell[0]
		    + celltranslation[1] * superCell[1]
		    + celltranslation[2] * superCell[2];
		  if (Length2(pos1 - pos[j->first]) < rCut2)
		    throw AsapError("Inconsistent neighbor list data: Did not find ")
		      << *i << " on normal nb list of atom " << j->first;
		}
#endif // CHECKCONSISTENCY
	    }
	}
    }

  // Get new neighbor lists for the modified atoms and for their new neighbors
  set<int> newneighbors;
  for (set<int>::const_iterator i = modified.begin(); i != modified.end(); ++i)
    {
      int l = cells->GetListAndTranslations(*i, nbList[*i]);
      l += cells->GetComplementaryListAndTranslations(*i, complNbList[*i]);
      if (l > maxLength)
	maxLength = l;
      for (nbiterat j = nbList[*i].begin(); j != nbList[*i].end(); ++j)
	if (!modified.count(j->first))
	  newneighbors.insert(j->first);
      for (nbiterat j = complNbList[*i].begin(); j != complNbList[*i].end();
	   ++j)
	if (!modified.count(j->first))
	  newneighbors.insert(j->first);
    }
  for (set<int>::const_iterator i = newneighbors.begin();
       i != newneighbors.end(); ++i)
    {
      int l = cells->GetListAndTranslations(*i, nbList[*i]);
      l += cells->GetComplementaryListAndTranslations(*i, complNbList[*i]);
      if (l > maxLength)
	maxLength = l;
    }
  affected.insert(newneighbors.begin(), newneighbors.end());
  // Now, we need to update the reference positions
  cells->UpdateReferencePositions(modified);

#if 0
  cerr << "Modified:";
  for (set<int>::const_iterator i = modified.begin();
       i != modified.end(); ++i)
    cerr << " " << *i;
  cerr << endl << "Affected:";
  for (set<int>::const_iterator i = affected.begin();
       i != affected.end(); ++i)
    cerr << " " << *i;
  cerr << endl;
    
#endif // 0
#ifdef CHKFULLLIST
  CheckFullListConsistency("End of RemakeLists");
#endif
}

int NeighborList::TestPartialUpdate(set<int> modified, PyObject *pyatoms)
{
  set<int> dummy;
  Atoms *atoms = GetAtoms();
  atoms->Begin(pyatoms);
  RemakeLists(modified, dummy);
  atoms->End();
  return dummy.size();
}

// Check consistency between nbList and complNbList.  Note that with
// Monte Carlo optimizations, consistency is not guaranteed for atoms
// separated by more than rCut.  If the distance is greater than
// rCut+2*drift the inconsistency is always safe, if lower than that
// but still greater than rCut is is most likely safe, but cannot be
// determined easily.

void NeighborList::CheckFullListConsistency(const string where, bool chkdst)
{
  RETURNIFASAPERROR;
  if (invalid || !fulllists)
    return;
  const vector<Vec> &pos = cells->GetWrappedPositions();
  const Vec *superCell = atoms->GetCell();
  double rc2 = rCut + 4*drift;
  double rCut2 = rCut * rCut;
  rc2 *= rc2;
  typedef vector< pair<int,translationsidx_t> >::const_iterator iterat;
  // Check that all info on nbList is found on complNbList
  for (int i = 0; i < nAtoms; i++)
    for (iterat j = nbList[i].begin(); j != nbList[i].end(); ++j)
      {
	// j->first is on the nblist of i.  Is i on complnblist of j->first?
	iterat k = complNbList[j->first].begin();
	iterat k_end = complNbList[j->first].end();
	while (k->first != i && k != k_end)
	  ++k;
	if (k->first != i)
	  {
	    // Is the inconsistency safe?
	    IVec celltranslation = translationTable[j->second];
	    Vec pos1 = pos[i] + celltranslation[0] * superCell[0]
	      + celltranslation[1] * superCell[1]
	      + celltranslation[2] * superCell[2];
	    if (Length2(pos1 - pos[j->first]) < rCut2)
	      throw AsapError("nbList[") << i << "] contains " << j->first
					 << ", but complNbList[" << j->first
					 << "] is missing " << i
					 << ". (" << where << ")";
	  }
	if (chkdst) {
	  // Are the atoms actually neighbors?
	  IVec celltranslation = translationTable[j->second];
	  Vec diff = pos[j->first] - pos[i] 
	    - celltranslation[0] * superCell[0]
	    - celltranslation[1] * superCell[1]
	    - celltranslation[2] * superCell[2];
	  if (Length2(diff) > rc2 + 1e-5)
	    throw AsapError("nbList[") << i << "] contains " << j->first
				       << ", but they are too far apart. ("
				       << where << ") "
				       << Length2(diff) << " > " << rc2
				       << "\npos1 = " << pos[i]
				       << "\npos2 = " << pos[j->first]
				       << "\ntranslation = " << celltranslation;
	}
      }
  // Check that all info on complNbList is found on nbList
  for (int i = 0; i < nAtoms; i++)
    for (iterat j = complNbList[i].begin(); j != complNbList[i].end(); ++j)
      {
	// j->first is on the complnblist of i.  Is i on nblist of j->first?
	iterat k = nbList[j->first].begin();
	iterat k_end = nbList[j->first].end();
	while (k->first != i && k != k_end)
	  ++k;
	if (k->first != i)
	  {
	    // Is the inconsistency safe?
	    IVec celltranslation = translationTable[j->second];
	    Vec diff = pos[j->first] - pos[i] 
	      - celltranslation[0] * superCell[0]
	      - celltranslation[1] * superCell[1]
	      - celltranslation[2] * superCell[2];
	    if (Length2(diff) < rCut2)
	      throw AsapError("complNbList[") << i << "] contains " << j->first
					      << ", but nbList[" << j->first
					      << "] is missing " << i
					      << ". (" << where << ")";
	  }
      }
}

void NeighborList::printlist(int n) const
{
  typedef vector< pair<int,translationsidx_t> >::const_iterator iterat;
  cerr << "nbList[" << n << "]";
  for (iterat i = nbList[n].begin(); i != nbList[n].end(); ++i)
    cerr << " " << i->first;
  cerr << endl;
  if (fulllists)
    {
      cerr << "complNbList[" << n << "]";
      for (iterat i = complNbList[n].begin(); i != complNbList[n].end(); ++i)
	cerr << " " << i->first;
      cerr << endl;
    }
} 

void NeighborList::print_info(int n)
{
  cerr << "NeighborList info on atom " << n << ":" << endl;
  cerr << "nbList:";
  for (int i = 0; i < nbList[n].size(); i++)
    cerr << "  " << nbList[n][i].first << " " << nbList[n][i].second;
  cerr << endl;
  if (fulllists)
    {
      cerr << "complNbList:";
      for (int i = 0; i < complNbList[n].size(); i++)
	cerr << "  " << complNbList[n][i].first << " "
	     << complNbList[n][i].second;
      cerr << endl;
    }
  cells->print_info(n);
}

long NeighborList::PrintMemory() const
{
  long n = 0;
  long ntot = 0;
  typedef vector< vector< pair<int, translationsidx_t> > >::const_iterator iter;
  for (iter i = nbList.begin(); i != nbList.end(); i++)
    {
      n += i->size();
      ntot += i->capacity();
    }
  if (fulllists)
    for (iter i = complNbList.begin(); i != complNbList.end(); i++)
      {
	n += i->size();
	ntot += i->capacity();
      }
  long mem = ntot * (sizeof(int) + sizeof(translationsidx_t));
  mem = (mem + 512*1024)/(1024*1024);
  long overhead = (ntot - n) * (sizeof(int) + sizeof(translationsidx_t));
  overhead = (overhead + 512*1024)/(1024*1024);
  char buffer[500];
  snprintf(buffer, 500,
	   "*MEM* NeighborList %ld MB.  [ overhead %ld MB, %.2e items, full=%d, sizeof(xlat_idx)=%ld ]",
	   mem, overhead, (double) n, (int) fulllists,
	   (long) sizeof(translationsidx_t));
  cerr << buffer << endl;
  mem += cells->PrintMemory();
  return mem;
}

double NeighborList::GetMaxStrainDisplacement()
{
  double disp2 = 0.0;
  const Vec *ss = atoms->GetCell();
  vector<Vec> strain(3);
  Vec factor;
  for (int i =0; i < 3; i++)
    {
      double length = sqrt(Length2(ss[i]));
      factor[i] = rCut / length;
      strain[i] = ss[i] - referenceSuperCell[i];
    }
  for (int i = -1; i <= 1; i++)
    for (int j = -1; j <= 1; j++)
      for (int k = -1; k <= 1; k++)
        {
          Vec v = i * factor[0] * strain[0] + j * factor[1] * strain[1] + k * factor[2] * strain[2];
          double l2 = Length2(v);
          if (l2 > disp2)
            disp2 = l2;
        }
  return sqrt(disp2);
}
