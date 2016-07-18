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
#include "RegularGridDecomposition.h"
#include "ParallelAtoms.h"
#include "AsapMPI.h"
#include "Exception.h"
#include <math.h>
#include <iostream>
#include <iomanip>
using std::cerr;
using std::endl;
using std::setw;

extern int verbose;

static const char offsetname[] = "_asap_rgd_offset";

RegularGridDecomposition::RegularGridDecomposition(const Vec superCell[3],
						   const bool periodic[3],
						   const int nCells[3], 
						   Communicator *comm)
{
  // Store information about the supercell and the CPU layout.
  nProcs = 1;
  for (int i = 0; i < 3; i++)
    {
      this->nCells[i] = nCells[i];
      nProcs *= nCells[i];
      originalPeriodicity[i] = periodic[i];
      this->superCell[i] = superCell[i];
    }
  nTotalCells[0] = 1;
  for (int i = 0; i < 3; i++)
    nTotalCells[i + 1] = nTotalCells[i] * nCells[i];
  this->comm = comm;
  assert(nProcs == comm->GetNumberOfProcessors());
  thisProcessor = comm->GetProcessorNumber();
  InitializeNeighborsTable();
#if 0
  if (thisProcessor == 0)
    {
      cerr << "RegularGridDecomposition neighbor table:" << endl;
      for (int i = 0; i < 27; i++)
        {
          cerr << setw(5) << i << ":";
          for (vector< pair<int, int> >::const_iterator j =
                 neighbors[i].begin(); j != neighbors[i].end(); ++j)
            cerr << " (" << j->first << "," << j->second << ")";
          cerr << endl;
        }
    }
#endif
  makeSendRecvLists();
}

void RegularGridDecomposition::whichProcessor(ParallelAtoms *atoms,
					      vector< vector<int> > &processor,
					      vector<int> &migrateAway)
{
  Vec minimum;
  Vec size;
  int nAtoms = atoms->GetNumberOfAtoms();
  vector<Vec> positions(nAtoms);
  const bool *pbc = atoms->GetBoundaryConditions();
  
  // Get the positions, and transform to scaled space.  Keep the
  // offset in a numpy array, so it can migrate with the atoms and be
  // reused when ghost export lists are calculated.
  npy_intp dims[2];
  dims[0] = nAtoms;
  dims[1] = 3;
  atoms->GetScaledPositions(positions);
  if (!atoms->AllFreeBoundaries())
    {
      PyObject *py_soffset;
      if (atoms->AllPeriodicBoundaries())
	py_soffset = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
      else
	py_soffset = PyArray_ZEROS(2, dims, NPY_DOUBLE, 0);
      if (py_soffset == NULL)
	throw AsapError("RegularGridDecomposition failed to create soffset");
      Vec *soffset = (Vec *) PyArray_DATA((PyArrayObject *) py_soffset);
      for (int i = 0; i < 3; i++)
	{
	  if (pbc[i])
	    for (int j = 0; j < nAtoms; j++)
	      {
		soffset[j][i] = floor(positions[j][i]);
		positions[j][i] -= soffset[j][i];
	      }
	}
      atoms->SetData(offsetname, py_soffset);
      Py_DECREF(py_soffset);
    }
  // Initialize the processor lists
  processor.resize(nProcs);
  for (int i = 0; i < nProcs; i++)
    processor[i].clear();
  migrateAway.clear();

  // Find the smallest box in scaled space that contains all the atoms
  smallestBox(positions, atoms->GetBoundaryConditions(),
	      nAtoms, minimum, size);
  // Store till later call to makeGhostExportList.
  this->minimum = minimum;
  this->size = size;
    
  // Find out which atoms belong to which processor:
  const Vec *x = &positions[0];
  for (int a = 0; a < nAtoms; a++)
    {
      int proc = 0;
      for (int i = 0; i < 3; i++)
	{
          int k = int(((*x)[i] - minimum[i]) / size[i] * nCells[i]);
          assert(k >= 0);
          if (k == nCells[i])
            k--;
          assert(k < nCells[i]);
          proc += nTotalCells[i] * k;
	}
      assert(proc >= 0);
      assert(proc < nProcs);
      x++;
      if (proc != thisProcessor)
        {
          processor[proc].push_back(a);
          migrateAway.push_back(a);
        }
    }
  if (verbose >= 2)
    {
      cerr << "RegularGridDecomposition (Node " << thisProcessor 
           << "): Sending to other processors: [";
      for (int i = 0; i < nProcs; i++)
        cerr << " " << processor[i].size();
      cerr << " ]" << endl;
    }
    
}

void RegularGridDecomposition::makeGhostExportLists(ParallelAtoms *atoms,
                                                    double rCut,
        vector< vector< pair<int, int> > > &ghostLists)
{
  Vec minimum = this->minimum;
  Vec size = this->size;
  int nAtoms = atoms->GetNumberOfAtoms();
  
  // Find the smallest box in scaled space that contains all atoms:
  Vec scaledRCut;
    
  for (int i = 0; i < 3; i++)
    {
      size[i] /= nCells[i];
    }

  const double *cellheights = atoms->GetCellHeights();
  for (int i = 0; i < 3; i++)
    minimum[i] += ((thisProcessor / nTotalCells[i]) % nCells[i]) * size[i];
  for (int i = 0; i < 3; i++)
    scaledRCut[i] = rCut / (cellheights[i] * size[i]);
  for (int i = 0; i < 3; i++)
    if (nCells[i] > 1 && scaledRCut[i] >= 0.5)
      throw AsapError("The system is too small for the number of processors along direction ") << i << ".";

  // Transform coordinates to scaled space:
  vector<Vec> positions(nAtoms);
  const bool *pbc = atoms->GetBoundaryConditions();
  atoms->GetScaledPositions(positions);
  // We must repeat the wrapping done when atoms were assigned to processors
  if (!atoms->AllFreeBoundaries())
    {
      PyObject *py_soffset = atoms->GetData(offsetname);
      Vec *soffset = (Vec *) PyArray_BYTES((PyArrayObject *) py_soffset);
      for (int i = 0; i < 3; i++)
	{
	  if (pbc[i])
	    for (int j = 0; j < nAtoms; j++)
	      positions[j][i] -= soffset[j][i];
	}
      Py_DECREF(py_soffset);
      atoms->DeleteData(offsetname);
    }

  if (verbose > 1)
    cerr << thisProcessor << ':' << minimum << " <> " << minimum + size
         << endl
         << nAtoms << " atoms on "
         << nCells[0] << " * " << nCells[1] << " * " << nCells[2] 
         << " processors" << endl;

  // Find out which atoms are ghost atoms on which processors:
  vector<Vec>::iterator p = positions.begin();
  for (int a = 0; a < nAtoms; a++)
    {
      int n = 1 + 3 + 3 * 3;
      int dn = 1;
      for (int i = 0; i < 3; i++, dn *= 3)
        if (nCells[i] > 1)
          {
            double s = ((*p)[i] - minimum[i]) / size[i];
            if (s < scaledRCut[i])
              n -= dn;
            else if (s > 1.0 - scaledRCut[i])
              n += dn;
          }
      typedef vector< pair<int, int> >::const_iterator VP;
      for (VP nb = neighbors[n].begin(); nb != neighbors[n].end(); ++nb)
        ghostLists[nb->first].push_back(pair<int, int>(a, nb->second));
      ++p;
    }
  assert(p == positions.end());
}

    
void RegularGridDecomposition::smallestBox(vector<Vec> &positions,
					   const bool* originalPeriodicity,
                                           int nAtoms, Vec &minimum, Vec &size)
{
  // Find the smallest box in scaled space that contains all atoms:
  for (int i = 0; i < 3; i++)
    {
      if (originalPeriodicity[i])
        {
          size[i] = 1.0;
          minimum[i] = 0.0;
#if 0
          double *p = (nAtoms > 0) ? &(positions[0][i]) : 0;
          for (int a = 0; a < nAtoms; a++, p += 3)
            {
              *p -= floor(*p);   // Can undoubtedly be done more efficient!
#if 0
              if (*p < 0.0 || *p > 1.0)
                cerr << "WARNING: position[" << a << "][" << i << "] = "
                     << setprecision(25) << setw(25) << *p << endl;
#endif
            }
#endif
        } 
      else
        {
          double min, max;
          if (nAtoms > 0)
            {
              double *p = &(positions[0][i]);
              min = *p;
              max = min;
              p += 3;
              for (int a = 1; a < nAtoms; a++, p += 3)
                if (*p > max)
                  max = *p;
                else if (*p < min)
                  min = *p;
            }
          else
            {
              min = 1e99;
              max = -1e99;
            }
          min = comm->Min(min);
          max = comm->Max(max);
          minimum[i] = min;
          size[i] = max - min;
          assert(max > -1e99 && min < 1e99);
          if (size[i] < 1e-7)
            {
              minimum[i] -= 0.5 * (1e-7 - size[i]);
              size[i] = 1e-7;
            }
        }
    }
    
  if (verbose > 1)
    cerr << "Node " << thisProcessor << ':' << minimum << " "
         << minimum + size << "   "
         << nAtoms << " atoms on "
         << nCells[0] << " * " << nCells[1] << " * " << nCells[2] 
         << " processors" << endl;
}

static inline int min(int x, int y) {return x < y ? x : y;}
static inline int max(int x, int y) {return x > y ? x : y;}

void RegularGridDecomposition::InitializeNeighborsTable()
{
  int cell[3];
  for (int i = 0; i < 3; i++)
      cell[i] = (thisProcessor / nTotalCells[i]) % nCells[i];
  int n = 0;
  int c[3];
  const bool *periodic = originalPeriodicity;
#if 0
  if (thisProcessor == 0)
    {
      cerr << "Periodic: " << periodic[0] << " " << periodic[1] << " "
	   <<periodic[2] << endl;
      cerr << "nCells: " << nCells[0] << " " << nCells[1] << " "
	   << nCells[2] << endl;
      cerr << "nTotalCells: " << nTotalCells[0] << " " << nTotalCells[1] << " "
	   << nTotalCells[2] << endl;
      cerr << "cell: " << cell[0] << " " << cell[1] << " "
	   << cell[2] << endl;
    }
#endif
  for (int j2 = -1; j2 <= 1; j2++)
  for (int j1 = -1; j1 <= 1; j1++)
  for (int j0 = -1; j0 <= 1; j0++, n++)
    for (c[2] = cell[2] + min(0, j2); c[2] <= cell[2] + max(0, j2); c[2]++)
    for (c[1] = cell[1] + min(0, j1); c[1] <= cell[1] + max(0, j1); c[1]++)
    for (c[0] = cell[0] + min(0, j0); c[0] <= cell[0] + max(0, j0); c[0]++)
      {
        int k[3];
        int i = 0;
        for (; i < 3; i++)
          {
            k[i] = c[i];
            if (periodic[i])
              {
                if (c[i] < 0)
                  k[i] = nCells[i] - 1;
                else if (c[i] >= nCells[i])
                  k[i] = 0;
              }
            else
              if (c[i] < 0 || c[i] >= nCells[i])
                break;
          }
        if (i < 3)
          continue;
        int p = k[0] + nTotalCells[1] * k[1] + nTotalCells[2] * k[2];
        if (p == thisProcessor)
          continue;
        int dn = 1;
        int nTranslation = 1 + 3 + 3 * 3;
        for (i = 0; i < 3; i++, dn *= 3)
          if (periodic[i])
          {
            if (c[i] != k[i] && nCells[i] == 1)
              break;
            nTranslation += (k[i] - c[i]) / nCells[i] * dn;
          }
        if (i < 3)
          continue;         
        neighbors[n].push_back(pair<int, int>(p, nTranslation));
      }
}

void RegularGridDecomposition::makeSendRecvLists()
{
  sendlist.clear();
  recvlist.clear();
  // Find the position of this processor in the grid
  int n = thisProcessor;
  int selfi = n % nCells[0];
  n = n / nCells[0];
  int selfj = n % nCells[1];
  int selfk = n / nCells[1];
  if (verbose)
    cerr << "I am processor " << thisProcessor << " at location "
	 << selfi << "," << selfj << "," << selfk << endl;
  // Loop over neighboring positions, +/- 1 in each direction if the grid is
  // big enough
  int imin = nCells[0] > 2 ? -1 : 0;
  int imax = nCells[0] > 1 ? 1 : 0;
  int jmin = nCells[1] > 2 ? -1 : 0;
  int jmax = nCells[1] > 1 ? 1 : 0;
  int kmin = nCells[2] > 2 ? -1 : 0;
  int kmax = nCells[2] > 1 ? 1 : 0;
  for (int is = imin; is <= imax; is++) {
    int ir = -is;
    if ((ir == -1) && (nCells[0] == 2)) ir = 1;
    for (int js = jmin; js <= jmax; js++) {
      int jr = -js;
      if ((jr == -1) && (nCells[1] == 2)) jr = 1;
      for (int ks = kmin; ks <= kmax; ks++) {
        int kr = -ks;
        if ((kr == -1) && (nCells[2] == 2)) kr = 1;
        // Get the processor numbers of the processors
        int node_s[3];
        int node_r[3];
        node_s[0] = selfi + is;
        node_s[1] = selfj + js;
        node_s[2] = selfk + ks;
        node_r[0] = selfi + ir;
        node_r[1] = selfj + jr;
        node_r[2] = selfk + kr;
        for (int i = 0; i < 3; i++) {
          if (node_s[i] < 0) node_s[i] += nCells[i];
          if (node_s[i] >= nCells[i]) node_s[i] -= nCells[i];
          if (node_r[i] < 0) node_r[i] += nCells[i];
          if (node_r[i] >= nCells[i]) node_r[i] -= nCells[i];
        }
        int sendto = nTotalCells[0] * node_s[0]
          + nTotalCells[1] * node_s[1] + nTotalCells[2] * node_s[2];
        int recvfrom = nTotalCells[0] * node_r[0]
          + nTotalCells[1] * node_r[1] + nTotalCells[2] * node_r[2];
        assert((sendto >= 0) && (recvfrom >= 0));
        assert((sendto < nProcs) && (recvfrom < nProcs));
        if ((is != 0) || (js != 0) || (ks != 0)) {
          sendlist.push_back(sendto);
          recvlist.push_back(recvfrom);
        } else {
	  if (verbose)
	    cerr << "DEBUG: sendto=" << sendto << " is=" << is << " js="
		 << js << " ks=" << ks << endl;
          assert(sendto == thisProcessor);
        }
      }
    }
  }
  if (verbose) {
    cerr << "Send list:    ";
    for (vector<int>::iterator i = sendlist.begin(); i != sendlist.end(); ++i)
      cerr << *i << " ";
    cerr << endl;
    cerr << "Receive list: ";
    for (vector<int>::iterator i = recvlist.begin(); i != recvlist.end(); ++i)
      cerr << *i << " ";
    cerr << endl;
  }
}

void RegularGridDecomposition::GetTranslations(ParallelAtoms *atoms,
					       vector<Vec> &translations) const
{
  const Vec *vectors = atoms->GetCell();
  int s = 0;
  int k[3];
  for (k[2] = -1; k[2] <= 1; k[2]++)
    for (k[1] = -1; k[1] <= 1; k[1]++)
      for (k[0] = -1; k[0] <= 1; k[0]++)
	{
	  Vec vec(0.0, 0.0, 0.0);
	  for (int i = 0; i < 3; i++)
	    vec += vectors[i] * k[i];
	  translations[s++] = vec;
	}
}
