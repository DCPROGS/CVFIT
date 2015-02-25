import os
import math
import numpy as np
from pylab import *

def plot_pooled(fit, axes=None, plotFit=False, legend=False,
    save_ASCII_name=None):

    axes.clear()
    axes.grid(True)

    axes.errorbar(fit.data.avX,
        fit.data.avY,
        yerr=fit.data.avS, fmt='o', ecolor='b', label='average')
    logplotX = np.log10(fit.data.avX)
    plotX = 10 ** np.linspace(np.floor(np.amin(logplotX) - 1),
        np.ceil(np.amax(logplotX)), 100)
    plotYg = fit.eq.equation(plotX,
        fit.eq.pars)
    axes.semilogx(plotX, plotYg, 'b-')

    if save_ASCII_name:
        name1 = save_ASCII_name[:-4] + '1.txt'
        name1.replace('.', '1.')
        fout = open(name1,'w')
        for i in range(len(fit.data.avX)):
            fout.write('{0:.6e}\t{1:.6e}\t{2:.6e}\n'.
                format(fit.data.avX[i], fit.data.avY[i], fit.data.avS[i]))
        fout.close()
        print 'file 1 saved'

    if plotFit:
        logplotX = np.log10(fit.data.X)
        plotX = 10 ** np.linspace(np.floor(np.amin(logplotX) - 1),
            np.ceil(np.amax(logplotX)), 100)
        plotYg = fit.eq.equation(plotX, fit.eq.pars)
        axes.semilogx(plotX, plotYg, 'b-')
        if save_ASCII_name:
            name2 = save_ASCII_name[:-4] + '2.txt'
            fout = open(name2,'w')
            for i in range(len(plotX)):
                fout.write('{0:.6e}\t{1:.6e}\n'.
                    format(plotX[i], plotYg[i]))
            fout.close()
            print 'file 2 saved'

    if legend:
        axes.legend(loc=2)






def plot(parent, axes=None, plotGuesses=False, plotFit=False, plotNorm=False, 
    pooled=False, legend=False): #, save_fig_name=None):

    axes.clear()
    axes.grid(True)

    if plotNorm:
        for session in parent.fits:
            axes.semilogx(session.data.X, 
                session.data.normY, 'o', label=session.data.title)
            logplotX = np.log10(session.data.X)
            plotX = 10 ** np.linspace(np.floor(np.amin(logplotX) - 1),
                np.ceil(np.amax(logplotX)), 100)
            plotYg = session.eq.equation(plotX, session.eq.normpars)
            axes.semilogx(plotX, plotYg, 'b-')

    elif pooled:
#            self.parent.canvas.axes.plot(self.parent.pooledfit.data.avX,
#                self.parent.pooledfit.data.avY, 'ro', label='average')
        axes.errorbar(parent.pooledfit.data.avX, 
            parent.pooledfit.data.avY, 
            yerr=parent.pooledfit.data.avS, fmt='o', ecolor='b', label='average')
        logplotX = np.log10(parent.pooledfit.data.avX)
        plotX = 10 ** np.linspace(np.floor(np.amin(logplotX) - 1),
            np.ceil(np.amax(logplotX)), 100)
        plotYg = parent.pooledfit.eq.equation(plotX,
            parent.pooledfit.eq.pars)
        axes.semilogx(plotX, plotYg, 'b-')

    else:
        for set in parent.data:
            if set.S.any() == 0:
                axes.semilogx(set.X, set.Y, 'o', label=set.title)
            else: 
                axes.errorbar(set.X, set.Y, yerr=set.S,
                    fmt='o', label=set.title)
                axes.set_xscale('log')

    if plotGuesses:
        for session in parent.fits:
            logplotX = np.log10(session.data.X)
            plotX = 10 ** np.linspace(np.floor(np.amin(logplotX) - 1),
                np.ceil(np.amax(logplotX)), 100)
            plotYg = session.eq.equation(plotX, session.eq.pars)
            axes.semilogx(plotX, plotYg, 'y-')

    if plotFit:
        for session in parent.fits:
            logplotX = np.log10(session.data.X)
            plotX = 10 ** np.linspace(np.floor(np.amin(logplotX) - 1),
                np.ceil(np.amax(logplotX)), 100)
            plotYg = session.eq.equation(plotX, session.eq.pars)
            axes.semilogx(plotX, plotYg, 'b-')

    if legend:
        axes.legend(loc=2)

def plot_hill_fit_result_single(filename, set, equation, 
    plotdata=True, plotfit=True, plotguess=False, plotaverage=False,
    save_fig=False, save_name=None):
    '''
    Plot the original data with guess and fitted curve in log graph.
    '''

    fig, ax = subplots()
    ax.set_xscale('log')
    ax.set_title(os.path.split(filename)[1][:-4] + ' ' + set.title)
    
    # Plot original data
    if plotdata:
        ax.plot(set.X, set.Y, 'ko', label='experimental data')

    # Plot average data.
    if plotaverage:
        ax.plot(set.avX, set.avY, 'ro', label='average')
        ax.errorbar(set.avX, set.avY, yerr=set.avS, fmt='none', ecolor='r')

    # Plot guess curve
    if plotguess:
        logplotX = np.log10(set.X)
        plotX = 10 ** np.linspace(np.floor(np.amin(logplotX) - 1),
                                  np.ceil(np.amax(logplotX) + 1), 100)
        plotYg = equation.equation(plotX, equation.guess)
        ax.plot(plotX, plotYg, 'y-', label='guess')
    # Plot fitted curve
    if plotfit:
        logplotX = np.log10(set.X)
        plotX = 10 ** np.linspace(np.floor(np.amin(logplotX) - 1),
                                  np.ceil(np.amax(logplotX) + 1), 100)
        plotYr = equation.equation(plotX, equation.pars)
        ax.plot(plotX, plotYr, 'b-', label='fit')

    legend(loc=2)
    
    if save_fig:
        fig.savefig(save_name)
        close(fig)
    else:
        show()

def plot_hill_fit_result_multiple(fname, fs, norm=False, save_fig=False, save_name=None):
    
    fig, ax = subplots()
    ax.set_xscale('log')
    type = ''
    if norm:
        type = ' (normalised)'
    ax.set_title(os.path.split(fname)[1][:-4]+type)
    
    # Generate the common x axis
    xmin = math.floor(math.log10(min([s.data.Xmin() for s in fs])))
    xmax = math.ceil(math.log10(max([s.data.Xmax() for s in fs])))
    plotx = 10 ** np.linspace(xmin, xmax, 100)
    
    for index, set in enumerate(fs):
        if norm:
            Y, pars = set.data.normY, set.eq.normpars
        else:
            Y, pars = set.data.Y, set.eq.pars
        # Plot original data and fitted curve
        ax.plot(set.data.X, Y, 'o', label=set.data.title)
        ax.plot(plotx, set.eq.equation(plotx, pars), 
            label=(set.data.title+' fit'))
                
    legend(loc=2)
    if save_fig:
        fig.savefig(save_name)
        close(fig)
    else:
        show()


