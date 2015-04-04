#! /usr/bin/python

__author__="remislp"
__date__ ="$29-Nov-2014 15:40:02$"

from cvfit import fitting
from cvfit import plots
from cvfit import data
from cvfit.fitting import SingleFitSession

if __name__ == "__main__":
    
    allsets, fname = fitting.load_data(example=True)
    print('File {0} loaded'.format(fname))
    print('{0:d} sets found.'.format(len(allsets)))
    # Ask which sets to use.
    sets = allsets
    
    sets = fitting.set_weights(sets)
    for i in range(len(sets)):
        print ('\nSet #{0:d}:'.format(i+1))
        print (sets[i])
    #settings = fitting.general_settings()
    
    eqname = fitting.choose_equation()
    print (eqname)
    if eqname == 'Hill' or eqname == 'Langmuir':
        from cvfit.equations import Hill as eqfit
         
    fitsessions = []
    for set in sets:
        equation = eqfit(eqname)
        fsession = SingleFitSession(set, equation)
        fsession.fit()
        fsession.calculate_errors()
        #plots.plot_hill_fit_result_single(fname, fsession.data, fsession.eq, plotguess=True)
        fitsessions.append(fsession)
    
    plots.plot_hill_fit_result_multiple(fname, fitsessions)
    pooldata = data.XYDataSet()
    for session in fitsessions:
        session.eq.normalise(session.data)
        pooldata.pool(session.data.X, session.data.normY, session.data.S)
        
    plots.plot_hill_fit_result_multiple(fname, fitsessions, norm=True)
    
    pooldata.weightmode = 1
    equation = eqfit(eqname)
    fsession = SingleFitSession(pooldata, equation)
    fsession.fit()
    fsession.calculate_errors()
    fsession.data.average_pooled()
    plots.plot_hill_fit_result_single(fname, fsession.data, fsession.eq,
        plotdata=False, plotaverage=True)
    