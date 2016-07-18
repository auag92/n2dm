from asap3 import *
from ase.lattice.cubic import FaceCenteredCubic
from ase.visualize.fieldplotter import *

clean = True
plotlog = "fieldplotter.log"

atoms = FaceCenteredCubic(size=(50,50,5), symbol="Ag")

def xy(a=atoms):
    r = atoms.get_positions()
    return r[:,0] + r[:,1]

plotter = FieldPlotter(atoms, xy, verbose=1)
plotter.set_dimensions((110,100))
plotter.set_data_range('plot')
plotter.set_output(PostScriptFile("fieldplot"))
plotter.set_output(GifFile("fieldplot"))
plotter.set_log(plotlog)
plotter.plot()

if clean:
    for name in (plotlog, "fieldplot0000.ps", "fieldplot0000.gif"):
        print "deleting", name
        os.unlink(name)


