#! /usr/bin/python

__author__="remislp"
__date__ ="$29-Nov-2014 15:40:02$"

from cvfit import fitting
from cvfit import plots
from cvfit.fitting import SingleFitSession

if __name__ == "__main__":
    
    allsets, fname = fitting.load_data(example=True)
    print('File {0} loaded'.format(fname))
    print('{0:d} sets found.'.format(len(allsets)))
    
    # Ask which sets to use.
    sets = allsets
    
    sets = fitting.set_weights(sets)
    for i in range(len(sets)):
        print '\nSet #{0:d}:'.format(i+1)
        print sets[i]
    
    #settings = fitting.general_settings()
    
    eqname = fitting.choose_equation()
    print eqname
    if eqname == 'Hill':
        from cvfit.hill import Hill as eqfit
        fixed = [True, False, False, False]
    if eqname == 'Langmuir':
        from cvfit.hill import Hill as eqfit
        fixed = [True, False, False, True]
        
    fitsessions = []
    for set in sets:
        
        equation = eqfit(eqname)
        equation.fixed = fixed
        
        fsession = SingleFitSession(set, equation)
        fsession.fit()
        fsession.calculate_errors()
        plots.plot_hill_guesses_result(fname, fsession.data, fsession.eq)
        fitsessions.append(fsession)
        
