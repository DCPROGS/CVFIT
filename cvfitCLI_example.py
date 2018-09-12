#! /usr/bin/python

__author__="remislp"
__date__ ="$29-Nov-2014 15:40:02$"

from cvfit import fitting
from cvfit import plots
from cvfit.fitting import SingleFitSession, MultipleFitSession

if __name__ == "__main__":
    
    # load data set
    sets, fname = fitting.load_data(example=True)
    print('File {0} loaded'.format(fname))
    print('{0:d} sets found.'.format(len(sets)))
    for i in range(len(sets)):
        print ('\nSet #{0:d}:'.format(i+1))
        print (sets[i])

    # load equation
    eqname = 'Hill'
    if eqname == 'Hill' or eqname == 'Langmuir':
        from cvfit.equations import Hill as eqfit
         
    # initiate fitting sessions, fit data and print results
    fitsessions = MultipleFitSession()
    for set in sets:
        equation = eqfit(eqname)
        fs = SingleFitSession(set, equation)
        fs.fit()
        fs.calculate_errors()
        print('\n*************************************************')
        print('\t' + fs.data.title + ' fit finished')
        print(fs.string_estimates())
        print(fs.string_liklimits())
        fitsessions.add(fs)
    print('\n*************************************************')
    print('\tAverage of all fits:')
    print(fitsessions.string_average_estimates())
    
    # plot fitted
    fplots = fitsessions.prepare_fplot('fit')
    fig = plots.cvfit_plot(sets, fig=None, 
        fplotsets=fplots, fplotline='b-',
        logX=True, logY=False, legend=True)
    
    # normalise to fitted maxima and replot
    for fs in fitsessions.list:
        fs.eq.normalise(fs.data)
    fplots = fitsessions.prepare_fplot('norm')
    plots.cvfit_plot(sets, fig=None, fplotsets=fplots, fplotline='b-',
        logX=True, logY=False, legend=True, norm=True)
    
    # pool data and fit pooled
    fitsessions.pool(norm=True)
    fitsessions.pooled.fit()
    fitsessions.pooled.calculate_errors()
    fitsessions.pooled.data.average_pooled()
    print('\n***********************\n**************************')
    print('\tNormalised and pooled data fit finished')
    print(fitsessions.pooled.string_estimates())
    print(fitsessions.pooled.string_liklimits())
    
    # plot pooled data fit
    fplots = fitsessions.prepare_fplot('pooled')
    plots.cvfit_plot([fitsessions.pooled.data], fig=None, fplotsets=fplots, 
        fplotline='b-', logX=True, logY=False, legend=False, pooled=True)
    