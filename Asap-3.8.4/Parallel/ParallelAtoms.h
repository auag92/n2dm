// -*- C++ -*-
// ParallelAtoms.h: Atoms access object for parallel simulations.
//
// Copyright (C) 2001-2011 Jakob Schiotz and Center for Individual
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
#include "Asap.h"
#include "Exception.h"
#include "Vec.h"
#include "SymTensor.h"
#include "Atoms.h"
#include <set>
#include <utility>
using std::set;
using std::pair;

namespace ASAPSPACE {

class Communicator;
// class RegularGridDecomposition;
class DomainDecomposition;

/// Extension of Atoms for parallel simulations.

/// The atoms in a parallel simulation are distributed on a regular
/// grid of processors.  The atoms on this processor are regular
/// atoms, atoms on nearby processors within the cutoff range of the
/// Potential are "ghost atoms".
class ParallelAtoms : public NormalAtoms
{
protected:
  virtual ~ParallelAtoms();

public:
  /// Create parallel atoms.
  ParallelAtoms(PyObject *py_atoms);

  // Reference counting
  friend void AsapAtoms_INCREF(Atoms *atoms);
  friend void AsapAtoms_DECREF(Atoms *atoms);

  /// Start accessing the Python atoms.
  ///
  /// Until End() is called, references are kept to the Python atoms
  /// and some of their attributes.
  virtual void Begin(PyObject *pyatoms, bool allow_reopen=false);

  /// Start accessing the Python atoms.
  ///
  /// Called with postmigrate=true *only* at the end of migration, where
  /// the atoms are closed and then reopened to update the number of
  /// atoms, but without triggering communication normally used to
  /// detect modified atoms.
  virtual void Begin(PyObject *pyatoms, bool allow_reopen, bool postmigrate);

  // End() is not redefined.


  /// Get total number of atoms in parallel simulation.
  ///
  /// In a serial simulation, same as GetNumberOfAtoms, in a parallel simulation
  /// it is the sum over these over all processors.
  virtual int GetTotalNumberOfAtoms() const {return nTotalAtoms;}

  /// Update flag and counters across processors.
  ///
  /// Called by a Potential with a flag indicating if the neighborlist
  /// should be updated.  In a serial simulation, just return the flag.
  ///
  /// In a parallel simulation, communicate across processors so the
  /// the flag passed as an argument is updated if just one processor
  /// thinks an update is necessary.  If the flag is true, a migration
  /// will also be triggered in a parallel simulation.
  virtual bool UpdateBeforeCalculation(bool flag, double range);

      /// Get a set of all elements present in the simulations.
  virtual void GetListOfElements(set<int> &elements) const;

  /// Distribute the atoms on the processors after creating them.

  /// Distribute should be called \em after the atoms have been
  /// created and any extra data in Python arrays (velocities etc) has
  /// been registered, but \em before the atoms are used for anything.
  void Distribute();

  /// Migrate the atoms between processors.
  void Migrate(bool distributing = false);

  /// Update data (positions, atomic numbers, ...) on ghost atoms.
  ///
  /// If copy_position is true, the ghost positions are copied into
  /// this object's positions array.  If false, this step is skipped
  /// as an optimization, this is done if the data will be copied
  /// anyway later.
  void UpdateGhostData();

  /// Create ghost atoms.
  void DecorateWithGhosts(double range);
  
  /// Get the Communicator object which does the actual communication.
  Communicator *GetCommunicator() {return mpi;}

  /// Return three integers specifying the CPU layout.
  virtual IVec GetNumberOfCells() const
  {return IVec(nCells[0], nCells[1], nCells[2]);}

  //  virtual void NormalizePositions(NeighborLocator *caller = 0);
  //  virtual void NormalizePositions(const set<int> &a, NeighborLocator *caller = 0)
  //  {throw AsapError("Internal error: Method ParallelAtoms::NormalizePositions(const set<int> &a) called");}
  
  const long *GetIdentities() const;
  
  /// Get the object which decides on which CPU the atoms belong.
  DomainDecomposition *GetDomainDecomposition() {return domainDecomp;}

  /// Get data for the ghost atoms.  Do we need these?
  //virtual const Vec *GetGhostPositions() const;
  //virtual const int *GetGhostAtomicNumbers() const;

  /// Communicate some kind of data to the ghost atoms.
  ///
  /// This is used if a potential needs to communicate some kind of
  /// date during the energy/force calculations.  If the data should
  /// be communicated at the beginning of the calculation, it should
  /// instead be placed in a ghost array.
  virtual void CommunicateData(double *address, int n = 1);

  /// Print memory usage
  virtual long PrintMemory() const;
  
  /// Receive and sum up data from ghost atoms (inverse communication).
  ///
  /// Used to sum up the forces.
  void CollectFromGhosts(vector<Vec> &data);
  void CollectFromGhosts(vector<SymTensor> &data);

private:
  /// Check that all atoms are present exactly once.
  void CheckIdentities();

  /// Take new boundary conditions into account.
  virtual void NewBoundaryConditions(const bool newperiodic[3]);
  
  /// Get the number of cells from the Python object.
  void extract_ncells(PyObject *pyatoms);

  /// Get a sorted list of array names.
  void get_array_names(PyObject *dict, vector<const char *> &names);

  /// Get a list of array objects.
  void get_arrays(PyObject *dict, vector<const char *> &names,
		  vector<PyArrayObject *> &arrays);

  /// Make a set of new arrays corresponding to a set of old arrays
  void make_new_arrays(vector<PyArrayObject *> &newarrays,
		       vector<PyArrayObject *> &oldarrays,
		       int size);

  /// Release references to arrays
  void release_arrays(vector<PyArrayObject *> &arrays, vector<const char *> &names,
		      int maxrefcount);

  /// Store arrays into a dictionary.
  void store_arrays(PyObject *dict, vector<const char *> &names,
		    vector<PyArrayObject *> &arrays);
  
  /// Store the max range where the ghosts are valid.  Negative: no ghosts yet.
  void set_ghost_range(double range);

  /// Read the max range where the ghosts are valid.  Negative: no ghosts yet.
  double get_ghost_range();

  /// Set the number of ghost atoms
  void set_number_of_ghosts(int nGhosts);
  
private:
  int nCells[3];                ///< number of cells (CPUs) in each direction
  int nProcessor;               ///< The number of this processor
  int nProcessors;              ///< The number of processors
  Communicator *mpi;            ///< This object does the actual communication.
  DomainDecomposition *domainDecomp; ///< Decides where atoms belong.
  int nTotalAtoms;              ///< Sum of the number of atoms on all nodes.
  int migrationCounter;         ///< Counts each time atoms migrate.
  vector<char> sendBuffer;      ///< Buffor for outbound communication
  vector<char> receiveBuffer;   ///< Buffer for inbound communication

  /// The ghost export list.
  ///
  /// This is a list of the ghost atoms on other processors which need
  /// to be updated from this one.
  ///
  /// ghosts is a list with one element for each processor.
  ///
  /// ghosts[i] is a list with one element for each atom on this
  /// processor that processor i needs to know about (it will be a
  /// ghost atom on processor i).
  ///
  /// ghosts[i][j] is a physical atom on this processor represented
  /// by a ghost atom on processor i.  ghosts[i][j] is a pair of
  /// two integers, the first is the index of the physical atom on
  /// this processor, the second is the index of the translation that
  /// should be applied to its position to give the position of the
  /// ghost.  The translation is non-zero if the other processor is a
  /// neighbor to this processor through the periodic boundaries.
  ///
  /// ghosts[p] is a vector of of ghost atoms that processor p needs
  /// to know about. A ghost is stored as a pair<int, int> = (n, t),
  /// so that the position of the ghost atom is equal to the position
  /// of atom number n translated with translation vector number t.
  ///
  /// In order to decide which atom should be a ghost atoms for which
  /// processors, each atom will get a index between 0 and 26 depending
  /// on where in the box it is positioned - or which boundary it is
  /// close to. In 2D the indices would run from 0 to 8 like this:
  ///
  /// \verbatim
  ///
  ///  ________________
  /// |6 |7         |8 |
  /// |__|__________|__|
  /// |3 |4         |5 |
  /// |  |          |  |
  /// |  |          |  |
  /// |  |          |  |
  /// |__|__________|__|
  /// |0 |1         |2 |
  /// |__|__________|__|
  ///
  /// \endverbatim
  ///
  /// The band arround the edge has a thickness of the cutoff plus
  /// twice the drift.
  vector< vector< pair<int, int> > > ghosts;

  int ghost_count;        ///< When where ghost positions last updated?
  bool decorated;         ///< Is the number of ghost atoms valid.
  bool justdecorated;     ///< First communication after DecorateWithGhosts?
  vector<int> num_ghosts_recv_from;
};

} // end namespace

