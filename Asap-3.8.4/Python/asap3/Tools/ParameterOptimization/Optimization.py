import sys
import time
import signal
import numpy as np

from exceptions import KeyboardInterrupt

from asap3.Tools.ParameterOptimization.Quantities import QuantityCalc
from asap3.Tools.ParameterOptimization.ScipyFmin import fmin, fmin_bfgs, fmin_powell
from asap3.Tools.Timing import Timer

class TimeoutException(Exception):
    pass

def timeout_handler(signum, frame):
    raise TimeoutException("Calculation exceeded the time limit!")

class ParameterOptimization():
    """Optimizes potential parameters based on given quantities.

    Input format of the potential parameters:

    {(elements): [parameters], ...}

    Input format of the quantities is:

    [[name, structure, elements, value, deviation], ...]

    The following quantities and structures are supported:

    Name                    Structures
    ----                    ----------
    cohesive_energy         fcc, hcp, bcc
    vacancy_energy          fcc, hcp, bcc
    heat_of_formation       l12, l10, B2
    phase_energy            fcc-hcp, fcc-bcc, bcc-hcp
    lattice_constant_a      fcc, hcp, l12, l10, B2
    lattice_constant_c      hcp, l10
    lattice_constant_ca     hcp, l10
    bulk_modulus            fcc, hcp, bcc, l12, B2
    elastic_constant_C11    fcc, hcp, bcc, l12, B2
    elastic_constant_C12    fcc, hcp, bcc, l12, B2
    elastic_constant_C13    hcp
    elastic_constant_C33    hcp
    elastic_constant_C44    fcc, hcp, bcc, l12, B2
    surface_energy          fcc100, fcc111, hcp0001
    surface_ratio           any combination of surface energies
                            on the form fcc100-fcc111
    impurity_energy         oct38-face, oct38-edge

    """

    def __init__(self, elements, calc_class, parameters, varparameters, quantities,
                 latticeconstants, calc_args=None, debug=False):
        self.calc_class = calc_class
        if isinstance(calc_args, tuple):
            self.calc_args = calc_args
        else:
            self.calc_args = ()

        self.debug = debug
        self.timer = Timer()

        # Interpret elements
        if not isinstance(elements, (tuple, list)):
            elements = [elements]
        #assert len(elements) < 3
        self.elements = list(elements)

        # Interpret parameters
        self.param_keys = []
        self.param_value = []
        self.param_variable = []
        for key in parameters.keys():
            if isinstance(key, tuple) and len(key) == 2:
                for symbol in key:
                    if symbol not in self.elements:
                        raise ValueError("Symbol in key must be in the element list")
                self.param_keys.append(key)
                self.param_value.append(parameters[key][:])
                self.param_variable.append(varparameters[key][:])
            else:
                raise KeyError("Parameter symbol key must be a tuple of lenght two.")

        # Interpret quantities
        self.quan_calc = QuantityCalc(latticeconstants, None, self.timer, False)
        self.quan_name = []
        self.quan_structure = []
        self.quan_elements = []
        self.quan_value = []
        self.quan_deviation = []
        for quan in quantities:
            if not hasattr(self.quan_calc, 'get_' + quan[0]):
                raise ValueError("The '%s' quantity cannot be calculated." % (quan[0],))
            self.quan_name.append(quan[0])
            self.quan_structure.append(quan[1])
            if isinstance(quan[2], (str, int)):
                self.quan_elements.append(tuple([quan[2]]))
            else:
                self.quan_elements.append(tuple(quan[2]))
            self.quan_value.append(quan[3])
            self.quan_deviation.append(quan[4])

    def valuelist_to_quantities(self, values, error=None):
        quantities = []
        for name, structure, elements, value in zip(self.quan_name,
                                                    self.quan_structure,
                                                    self.quan_elements,
                                                    values):
            key = (name, structure, elements)
            if error != None and error >= 1e9 - 2:
                value = 0.0
            quantities.append([name, structure, elements, value])

        return quantities

    def quantities_to_valuelist(self, quantities):
        values = []
        for i, quan in enumerate(quantities):
            assert self.quan_name[i] == quan[0]
            values.append(quan[3])
        return values

    def parameters_list2dict(self, paramlist, copyall=False):
        """Copies a list of parameters (used in optimization) to a dictionary
        (used in potentials). Copies by default only variable parameters - set
        copyall=True to copy all parameters."""
        k = 0
        paramdict = {}
        for i, key in enumerate(self.param_keys):
            paramdict[key] = self.param_value[i][:]
            for j, bool in enumerate(self.param_variable[i]):
                if bool or copyall:
                    paramdict[key][j] = paramlist[k]
                    k += 1
        return paramdict

    def parameters_dict2list(self, paramdict=None, copyall=False):
        """Reverse of parameters_list2dict."""
        paramlist = []
        for i, key in enumerate(self.param_keys):
            for j, bool in enumerate(self.param_variable[i]):
                if bool or copyall:
                    if paramdict == None:
                        paramlist.append(self.param_value[i][j])
                    else:
                        paramlist.append(paramdict[key][j])
        return paramlist

    def fit(self, parameters=None, xtol=1e-2, ftol=5e-2, delta=5e-3,
            log=None, err=None, rawout=False):
        # Open and print log
        if log is not None:
            if isinstance(log, str):
                self.log = open(log, 'a')
            else:
                # Hopefully an open file.
                self.log = log
            #assert hasattr(self.log, 'write')
            
            output = "Time      Error Function    "
            format = "%02i:%02i:%02i  %-18.8e"
            T = time.localtime()
            strargs = [T[3], T[4], T[5], 0.0]
            for i, key in enumerate(self.param_keys):
                for j, bool in enumerate(self.param_variable[i]):
                    if bool:
                        output += "P%i.%-8i" % (i,j)
                        format += "%-11.6f"
                        strargs.append(0.0)
            for i, value in enumerate(self.quan_value):
                output += "Q%-9s" % (i,)
                format += "%-10.3f"
                strargs.append(value)
            T = time.localtime()
            self.log.write("\n" + output + "\n")
            self.log.write(len(output) * "-" + "\n")
            self.log.write(format % tuple(strargs) + "\n")
            self.log.flush()
            self.log_format = format
        else:
            self.log = None

        # Open error file/pipe
        if err is not None:
            if isinstance(err, str):
                self.err = open(err, 'w')
            else:
                # Hopefully an open file.
                self.err = err
            #assert hasattr(self.err, 'write')
        else:
            self.err = None

        # Copy the variable parameters to a list - if its is a list assume 
        # only varible parameters
        if not isinstance(parameters, list):
            parameters = self.parameters_dict2list(parameters)

        # Start optimization
        xopt, fopt, direc, iter, calls, flag = fmin_powell(self.error_function, parameters,
                                                           xtol=xtol, ftol=ftol,
                                                           full_output=True, disp=False)

        # Print log after optimization
        if self.log is not None:
            if flag == 0:
                self.log.write(("\n\nFitting terminated successfully.\n"))
            else:
                self.log.write(("\n\nFitting not converged - " +
                                "too many iterations or function calls.\n"))
            self.log.write(("   Error function value: %f\n" % (fopt,) +
                            "   Iterations: %d\n" % (iter,) +
                            "   Error function evaluations: %d\n\n\n" % (calls,)))
            #self.timer.write(self.log)
        else:
            if flag > 0:
                print ("Fitting not converged - " +
                       "too many iterations or function calls!")

        # Replace the old parameters with the optimized ones
        parameters = self.parameters_list2dict(xopt)

        # Get last calculated values of quantities
        (error, values) = self.calc_values(parameters)
        quantities = self.valuelist_to_quantities(values, error)

        if rawout:
            return (error, list(xopt), values)
        else:
            return (error, parameters, quantities)

    def calc_values(self, parameters):
        # Set calculator and calculate values
        self.timer.start('calc_quantity_values')
        try:
            calc = self.calc_class(self.elements[:], parameters, *self.calc_args)
            self.quan_calc.set_calculator(calc)
        except:
            T = time.localtime()
            self.err.write("Calculation Error..." +
                           "\n   Time: %02i-%02i-%i" % (T[2], T[1], T[0]) +
                           " %02i:%02i:%02i" % (T[3], T[4], T[5]) +
                           "\n   Type: %s" % (sys.exc_info()[0],) + 
                           "\n   Message: %s" % (sys.exc_info()[1],) + 
                           "\n   Parameters: %s" % (parameters,) + "\n\n")
            self.err.flush()
            values = [0.0] * len(self.quan_name)
            error = 1e100 - 1
        else:
            error = 0.0
            values = []
            for name, structure, elements, fitvalue, fitdev in zip(self.quan_name,
                                                                   self.quan_structure,
                                                                   self.quan_elements,
                                                                   self.quan_value,
                                                                   self.quan_deviation):
                old_handler = signal.signal(signal.SIGALRM, timeout_handler)
                signal.alarm(2)
                # Triger an alarm in 2 seconds. Simpel statistics show that the
                # mean calculation time for a quantity lies around 0.5 seconds,
                # with a maximun registred around 1 second.

                try:
                    try:
                        value = self.quan_calc.get_quantity(name, structure, elements)
                        if fitvalue == 0:
                            error += (value - fitvalue)**2 / (fitdev)**2
                        else:
                            error += (value - fitvalue)**2 / (fitvalue * fitdev)**2
                    except KeyboardInterrupt:
                        raise
                    except: #(TimeoutException, RuntimeError, AsapError):
                        T = time.localtime()
                        self.err.write("Calculation Error..." +
                                       "\n   Time: %02i-%02i-%i" % (T[2], T[1], T[0]) +
                                       " %02i:%02i:%02i" % (T[3], T[4], T[5]) +
                                       "\n   Type: %s" % (sys.exc_info()[0],) + 
                                       "\n   Message: %s" % (sys.exc_info()[1],) + 
                                       "\n\n   Value: %s" % (name,) +
                                       "\n   Structure: %s" % (structure,) +
                                       "\n   Elements: %s" % (elements,) +
                                       "\n   Parameters: %s" % (parameters,) +
                                       "\n   Values: %s" % (values,) + "\n\n")
                        self.err.flush()
                        values = [0.0] * len(self.quan_name)
                        error = 1e100 - 1
                        break
                    else:
                        values.append(value)
                except TimeoutException:
                    pass
                finally:
                    signal.alarm(0)
                    signal.signal(signal.SIGALRM, old_handler)

        self.timer.stop('calc_quantity_values')
        return (error, values)

    def error_function(self, newparameters):
        # Replace the old parameters with the new ones
        parameters = self.parameters_list2dict(newparameters)

        # Calculate values and errorfunction
        (errfunc, values) = self.calc_values(parameters)

        # Print log
        if self.log is not None:
            #if errfunc > 1e9 - 1:
            #    errfunc = 1e9 - 1
            T = time.localtime()
            strargs = [T[3], T[4], T[5], errfunc]
            strargs += newparameters.tolist()
            strargs += values
            self.log.write(self.log_format % tuple(strargs) + "\n")
            self.log.flush()

        return errfunc

    def write_result(self, fitpar, optpar, optval, out=sys.stdout):
        if isinstance(fitpar, list):
            fitpar = self.parameters_list2dict(fitpar)
        if isinstance(optpar, list):
            optpar = self.parameters_list2dict(optpar)
        if isinstance(optval[0], (tuple, list)):
            optval = self.quantities_to_valuelist(optval)

        for i, key in enumerate(self.param_keys):
            out.write("  %s%s parameters " % (key[0], key[1]) +
                      "(Initial, Optimal, Variable):\n")
            for pi, po, v in zip(fitpar[key], optpar[key],
                                 self.param_variable[i]):
                out.write("  %12.9f  %12.9f  %s\n" % (pi, po, v))
        out.write("\n  Fitting values (Quantity, Structure: " +
                  "Target, Optimal, Weight, Deviation):\n")
        devsum = 0
        for i, name in enumerate(self.quan_name):
            elements = self.quan_elements[i]
            structure = self.quan_structure[i]
            value = self.quan_value[i]
            weight = self.quan_deviation[i]

            if len(elements) > 1:
                structure = "%s%s-%s" % (elements + (structure,))
            else:
                structure = "%s-%s" % (elements + (structure,))
            name = "%s, %s:" % (name, structure)
            if value == 0:
                dev = np.abs(optval[i]) * 100
            else:
                dev = np.abs(optval[i] / value - 1) * 100
            devsum += dev
            out.write("  %-45s  %7.3f  " % (name, value) +
                      "%10.6f  %5.3f  " % (optval[i], weight) +
                      "%7.2f%%\n" % (dev,))
        out.write("%84.2f%%\n\n" % (devsum / len(self.quan_name),))

if __name__ == '__main__':
    from asap3.Tools.ParameterOptimization.EMT import EMT2011Fit, EMTStdParameters

    # Parameter names: ['eta2', 'lambda', 'kappa', 'E0', 'V0', 'S0', 'n0']
    initparam = {('Pt','Pt'): EMTStdParameters('Pt')}
    varparam = {('Pt','Pt'): [True, True, True, True, True, True, False]}

    quantities = [('lattice_constant_a', 'fcc', 'Pt', 3.92,      0.001),
                  #('lattice_constant_a', 'hcp', 'Pt', 2.77,      0.05),
                  #('lattice_ratio_ca', 'hcp', 'Pt', 4.78 / 2.77, 0.05),
                  ('bulk_modulus', 'fcc', 'Pt', 278.3,           0.01),
                  #('elastic_anisotropy', 'fcc', 'Pt', 1.594,   0.01),
                  #('elastic_constant_C11', 'fcc', 'Pt', 346.7,   0.1),
                  #('elastic_constant_C12', 'fcc', 'Pt', 250.7,  0.1),
                  #('elastic_constant_C44', 'fcc', 'Pt', 76.5,    0.01),
                  ('cohesive_energy', 'fcc', 'Pt', 5.84,         0.001),
                  #('phase_energy', 'fcc-hcp', 'Pt', -0.05,       0.02),
                  #('surface_energy', 'fcc111', 'Pt', 0.631,      0.02),
                  #('surface_ratio', 'fcc111-fcc100', 'Pt', 0.631 / 0.892, 0.01),
                  ]

    latticeconstants = [['fcc', 'Pt', 3.92],
                        ['hcp', 'Pt', (2.77, 4.78)]]

    opt = ParameterOptimization(['Pt'], EMT2011Fit, initparam, varparam,
                                quantities, latticeconstants, ('kappa',), False)
    (error, optpar, optval) = opt.fit(xtol=0.1, ftol=10.0, delta=0.01,
                                      log=sys.stdout, err='Optimization.err')
    print 'Optimization ended with error function %.6e' % (error,)
    opt.write_result(initparam, optpar, optval)

