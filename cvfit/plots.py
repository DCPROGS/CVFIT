import os
import numpy as np
from pylab import *

def plot_hill_guesses_result(filename, set, equation):
    '''
    Plot the original data with guess and fitted curve in log graph.
    '''

    fig, ax = subplots()
    ax.set_xscale('log')
    ax.set_title(
    os.path.split(filename)[1][:-4] + ' ' + set.title)

    # Plot original data
    pl1, = ax.plot(set.X, set.Y, 'o', label='experimental data')
    
    logplotX = np.log10(set.X)
    plotX = 10 ** np.linspace(np.floor(np.amin(logplotX) - 1),
                              np.ceil(np.amax(logplotX) + 1), 100)
    #xmin = plotX[0].copy()
    #xmax = plotX[-1].copy()

    # Plot guess curve
    plotYg = equation.equation(plotX, equation.guess)
    pl2, = ax.plot(plotX, plotYg, 'r-', label='guess')
    # Plot fitted curve
    plotYr = equation.equation(plotX, equation.pars)
    pl3, = ax.plot(plotX, plotYr, 'b-', label='fit')
    
    legend(handles=[pl1, pl2, pl3], loc='lower right')
    show()


    fig.savefig(filename[:-4] + '_' + set.title + '.png')
    plt.close(fig)

