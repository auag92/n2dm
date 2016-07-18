# Copyright (C) 2003  CAMP
# Please see the accompanying LICENSE file for further information.

import os
import sys

import numpy as npy

from asap3 import _asap

MASTER = 0
debug = True

def is_contiguous(array, dtype=None):
    """Check for contiguity and type."""
    if dtype is None:
        return array.flags.contiguous
    else:
        return array.flags.contiguous and array.dtype == dtype

# Serial communicator
class SerialCommunicator:
    size = 1
    rank = 0
    def sum(self, array, root=-1):
        if isinstance(array, (float, complex)):
            return array

    def scatter(self, s, r, root):
        r[:] = s

    def max(self, value, root=-1):
        return value

    def broadcast(self, buf, root):
        pass

    def send(self, buff, root, tag=123, block=True):
        pass

    def barrier(self):
        pass

    def gather(self, a, root, b):
        b[:] = a

    def new_communicator(self, ranks):
        return self

    def cart_create(self, dimx, dimy, dimz, periodic):
        return self

class DummyCommunicator(SerialCommunicator):
    size = 1
    size = 0
    def new_communicator(self, ranks):
        new_comm = DummyCommunicator()
        new_comm.size = len(ranks)
        return new_comm


# serial_comm = SerialCommunicator()

try:
    world = _asap.Communicator()
except AttributeError:
    world = SerialCommunicator()

size = world.size
rank = world.rank
parallel = (size > 1)

if debug and parallel:
    class _Communicator(_asap.Communicator):
        def new_communicator(self, ranks):
            assert is_contiguous(ranks, int)
            sranks = npy.sort(ranks)
            # Are all ranks in range?
            assert 0 <= sranks[0] and sranks[-1] < self.size
            # No duplicates:
            for i in range(len(sranks) - 1):
                assert sranks[i] != sranks[i + 1]
            comm = _asap.Communicator.new_communicator(self, ranks)
            if comm is None:
                # This cpu is not in the new communicator:
                return None
            else:
                # Oops, this will be the original _asap.Communicator class.
                return comm

        def broadcast_string(self, string, root):
            if rank == root:
                assert isinstance(string, str)
                n = npy.array(len(string), int)
            else:
                n = npy.zeros(1, int)
            self.broadcast(n, root)
            if rank == root:
                string = npy.fromstring(string, npy.int8)
            else:
                string = npy.zeros(n, npy.int8)
            self.broadcast(string, root)
            return string.tostring()

    world = _Communicator()
    assert size == world.size and rank == world.rank

def broadcast_string(string=None, root=MASTER, comm=world):
    if rank == root:
        assert isinstance(string, str)
        n = npy.array(len(string), int)
    else:
        assert string is None
        n = npy.zeros(1, int)
    comm.broadcast(n, root)
    if rank == root:
        string = npy.fromstring(string, npy.int8)
    else:
        string = npy.zeros(n, npy.int8)
    comm.broadcast(string, root)
    return string.tostring()


def all_gather_array(comm, a): #???
    # Gather array into flat array
    shape = (comm.size,) + npy.shape(a)
    all = npy.zeros(shape)
    comm.all_gather(a, all)
    return all.ravel()


# def run(iterators):
#     """Run through list of iterators one step at a time."""
#     if not isinstance(iterators, list):
#         # It's a single iterator - empty it:
#         for i in iterators:
#             pass
#         return

#     if len(iterators) == 0:
#         return

#     while True:
#         try:
#             results = [iter.next() for iter in iterators]
#         except StopIteration:
#             return results
