# Copyright (C) 2014 Simon H. Brodersen and Center for Individual
# Nanoparticle Functionality, Department of Physics, Technical
# University of Denmark.  Email: schiotz@fysik.dtu.dk
#
# This file is part of Asap version 3.
#
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# version 3 as published by the Free Software Foundation.  Permission
# to use other versions of the GNU Lesser General Public License may
# granted by Jakob Schiotz or the head of department of the
# Department of Physics, Technical University of Denmark, as
# described in section 14 of the GNU General Public License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# and the GNU Lesser Public License along with this program.  If not,
# see <http://www.gnu.org/licenses/>.

import numpy as np

from asap3 import FullNeighborList, CoordinationNumbers
from asap3.analysis.localstructure import GuessLatticeConstant

def GetDiameter(atoms, N=500, all=False):
    # Based on the spiral points method by E. Saff and A. Kuijlaars,
    # The Mathematical Intelligencer (1997)
    h = 2.0*np.arange(N)/(N - 1) - 1
    theta = np.arccos(h)
    phi = np.zeros(N)
    for i in range(1, N - 1):
        phi[i] = (phi[i - 1] + 3.6/np.sqrt(N*(1 - h[i]**2)))# % 2*np.pi

    directions = np.empty((N,3))
    directions[:,0] = np.sin(theta) * np.cos(phi)
    directions[:,1] = np.sin(theta) * np.sin(phi)
    directions[:,2] = np.cos(theta)

    # Project the positions onto the grid and return the mean diameter
    p = atoms.get_positions()
    r = p - p.mean(axis=0)

    diameters = np.zeros(N)
    for i, d in enumerate(directions):
        d = np.dot(r, d)
        diameters[i] = d.max() - d.min()

    if all:
        return diameters
    else:
        return diameters.mean()

def GetLayerNumbers(atoms, rNearest=None, maxlayers=-1):
    if not rNearest:
        rNearest = GuessLatticeConstant(atoms) * 0.707
    c_nb = rNearest * 1.17 #Nearest neighbor cutoff
    c_sh = rNearest * 1.50 #Surface atoms cutoff 

    nb = FullNeighborList(c_nb, atoms)
    cn = CoordinationNumbers(atoms, c_nb)

    rad = GetDiameter(atoms) / 2.0
    pos = atoms.get_positions()
    cop = pos.mean(axis=0)
    r = np.sqrt(np.sum((pos - cop)**2, axis=1))

    tags = ((cn < 11) & (r > (rad - c_sh))).astype(int)
    left = len(atoms) - tags.sum()

    n = 1
    while left:
        for i in np.arange(len(tags))[tags == n]:
            nbi = nb[i]
            ids = nbi[tags[nbi] == 0]
            tags[ids] = n + 1
            left -= len(ids)
        n += 1
        if n > maxlayers and maxlayers > 0:
            tags[tags == 0] = n
            break

    return tags

