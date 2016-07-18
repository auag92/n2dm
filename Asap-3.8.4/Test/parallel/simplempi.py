import asap3
from asap3.mpi import world
from asap3.testtools import ReportTest
import numpy as np

if world.size == 1:
    raise RuntimeError("Cannot test MPI communication on a single processor.")

if world.rank == 0:
    asap3.print_version(1)

world.barrier()

outdata = np.array([world.rank * 2], int)
indata = np.zeros(1, int)

to = world.rank + 1
if to >= world.size:
    to = 0
fr = world.rank - 1
if fr < 0:
    fr += world.size

print "%d: to=%d from=%d" % (world.rank, to, fr)

rq = world.send(outdata, to, block=False)
world.receive(indata, fr)
world.wait(rq)

ReportTest("Data received by %s from %s" % (world.rank, fr), indata[0], 2*fr, 0)


outdata = world.rank * np.arange(10000000)
indata = np.zeros(10000000, int)

to = world.rank + 1
if to >= world.size:
    to = 0
fr = world.rank - 1
if fr < 0:
    fr += world.size

print "%d: to=%d from=%d" % (world.rank, to, fr)

rq = world.send(outdata, to, block=False)
world.receive(indata, fr)
#world.wait(rq)
print "Request completed:", rq.test()
rq.wait()

assert (indata == fr * np.arange(10000000)).all()

s = world.sum(2.5 * world.rank)
ReportTest("Float sum", s, 2.5 * 0.5 * world.size * (world.size - 1), 1e-9)

s = world.sum(world.rank)
ReportTest("Integer sum", s, world.size * (world.size - 1) / 2, 0)

s = world.max(1.0 * world.rank)
ReportTest("Float max", s, world.size - 1, 0)

s = world.max(world.rank)
ReportTest("Int max", s, world.size - 1, 0)

print "Checking reduces on arrays..."
d0 = np.arange(10) + world.rank
f0 = d0 + 0.0

d = d0.copy()
f = f0.copy()
world.max(d)
assert (d == np.arange(10) + world.size - 1).all()
world.max(f)
assert (f == np.arange(10) + world.size - 1.0).all()

d = d0.copy()
f = f0.copy()
world.min(d)
assert (d == np.arange(10)).all()
world.min(f)
assert (f == np.arange(10)).all()

print "  ... passed."

print "Checking new communicator"
newcomm = world.new_communicator(np.arange(2))
if newcomm:
    print "I am", newcomm.rank
    s = newcomm.sum(world.rank+1)
    ReportTest("Integer sum", s, 3, 0)

world.barrier()
ReportTest.Summary()
