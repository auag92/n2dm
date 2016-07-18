This folder contains:

AdsCalc.py: The calculator to be used instead of say EMT to account for the gas

AdsorptionParameter.py: Contains gas parameters for different materials

atommontecarlodata.py: The data class for the amc simulations

atomsmontecarlo.py: Starts and atom simulation with or without gas.

endatomsimulation.py: finishes an atom monte carlo simulation and writes result to log.

langmuirExpression.py: returns the gas coverages and energies(the only one to use directly)

Logger.py: Logging module.

queuemultipleamc.py: Queues amc simulation and resume scripts that depend on the preceeding.

resizecluster.py: corrects the n.o. atoms before amc simulation.

resume_amc_gas.py: resumes a prematurly stopped amc simulation.

surface_end.py: ends a surface simulation by filtering.

surface_start: submits a n.o. surface simulations.
