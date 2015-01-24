import os
import math
import numpy as np
from pylab import *
import prettyplotlib as ppl

def plot_hill_fit_result_single(filename, set, equation, 
    plotdata=True, plotfit=True, plotguess=False, plotaverage=False,
    save_fig=False, save_name=None):
    '''
    Plot the original data with guess and fitted curve in log graph.
    '''

    fig, ax = subplots()
    ax.set_xscale('log')
    ax.set_title(os.path.split(filename)[1][:-4] + ' ' + set.title)
    handles = []
    
    # Plot original data
    if plotdata:
        pl1, = ax.plot(set.X, set.Y, 'ko', label='experimental data')
        handles.append(pl1)
    # Plot average data.
    if plotaverage:
        pl4, = ax.plot(set.avX, set.avY, 'ro', label='average')
        ax.errorbar(set.avX, set.avY, yerr=set.avS, fmt='none', ecolor='r')
        handles.append(pl4)
    # Plot guess curve
    if plotguess:
        logplotX = np.log10(set.X)
        plotX = 10 ** np.linspace(np.floor(np.amin(logplotX) - 1),
                                  np.ceil(np.amax(logplotX) + 1), 100)
        plotYg = equation.equation(plotX, equation.guess)
        pl2, = ax.plot(plotX, plotYg, 'y-', label='guess')
        handles.append(pl2)
    # Plot fitted curve
    if plotfit:
        logplotX = np.log10(set.X)
        plotX = 10 ** np.linspace(np.floor(np.amin(logplotX) - 1),
                                  np.ceil(np.amax(logplotX) + 1), 100)
        plotYr = equation.equation(plotX, equation.pars)
        pl3, = ax.plot(plotX, plotYr, 'b-', label='fit')
        handles.append(pl3)
    
    legend(handles=handles, loc='lower right')
    
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
    
    handles = []
    for index, set in enumerate(fs):
        if norm:
            Y = set.data.normY
            pars = set.eq.normpars
        else:
            Y = set.data.Y
            pars = set.eq.pars
        # Plot original data and fitted curve
        pl1, = ax.plot(set.data.X, Y,
                 'o', color=ppl.colors.set2[index], label=set.data.title)
        pl2, = ax.plot(plotx, set.eq.equation(plotx, pars),
                 color=ppl.colors.set2[index], label=(set.data.title+' fit'))
        handles.extend([pl1, pl2])
    legend(handles=handles, loc='upper left')
    
    if save_fig:
        fig.savefig(save_name)
        close(fig)
    else:
        show()
