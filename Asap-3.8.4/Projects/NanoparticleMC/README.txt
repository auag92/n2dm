This folder contains script files for two-level monte carlo
simulations of nanoparticles.

LIST OF SCRIPTS IN ORDER OF USAGE.

All scripts can be run with -h to print a usage message.


surface_start.py:

Runs the surface monte carlo part.  Jobs are submitted to
Niflheim by the script.

Note on temperatures: temp is the fake simulation temperature, it
should be chosen high to sample phase space appropriately, e.g. 4000K.
tgas is the temperature used to calculate the gas adsorption.  It is a
real physical temperature, and should be set to the temperature at
which the final result is desired, e.g. 300K.


surface_end.py:

Harvests the results from the surface monte carlo runs, and collects
the lowest-energy configurations in a single .smc file.


atomsmontecarlo.py:

Runs the atoms-monte-carlo simulations.  The temperature is a
fictitious temperature that is higher than the physical temperature
but lower than the one used in the surface monte carlo.  Eg 1200K.
tgas is the real temperature, e.g. 300K


endatomssimulation.py:

Harvests the results from the atoms-monte-carlo simulations, and
converts to the real temperature (300K in this example).


SCRIPTS THAT HAVE NOT BEEN UPDATED:

queuemultipleamc.py
resume_amc_gas.py


TO DO:

Move parameters into a single parameter file, to ensure consistency
between simulations.



