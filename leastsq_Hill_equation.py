import sys
import os
import codecs


import numpy as np
from scipy import optimize, stats
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import markdown
import prettyplotlib as ppl


from toolkit import Report, check_input
from hillequation import HillEquation, hill_equation, hill_equation_residuals, \
    fitting_curve_hill_equation


class Cell(object):

    '''
    A class that store the inportant information for each cell
    '''

    def __init__(self, concentration, response):
        '''
        Make the class by specifying the concentration and response.
        Original data is then sorted according to the concentration.
        '''

        # Keep track of the original data
        self.originalconcentration = concentration
        self.originalresponse = response
        # Sort the data according to concentration
        result = np.vstack((concentration, response))
        result.sort()
        self.concentration = result[0]
        self.response = result[1]
        self.size = len(result[0])

print 'Please type in the loaction of the csv file'
filename = raw_input('filename:').strip()
while (os.path.exists(filename) is False) or (filename[-4:] != '.csv'):
    if os.path.exists(filename) is False:
        print filename, 'not found. Please type in again.'
    elif filename[-4:] != '.csv':
        print filename, 'is not a CSV file. Please type in again.'
    filename = raw_input('filename:').strip()

report = Report(filename[:-4])

rawresult = np.genfromtxt(filename, delimiter=',')
setnum = rawresult.shape[1] / 2
print str(setnum), 'sets found.'
report.title('Original data:', 1)
report.paragraph('Number of cells been examined: ' + str(setnum))

print 'Please type in the mode number:'
print 'Mode 1: Do different fitting for different data sets.(Default)'
print 'Mode 2: Do one fitting for all the data sets.'
rearange = check_input('Mode number:', ['1', '2'], 1)

if rearange == 2:
    result = rawresult[np.isfinite(rawresult)]
    celllist = [Cell(result[::2], result[1::2])]
elif rearange == 1:
    celllist = []
    for index in range(setnum):
        real = np.isfinite(rawresult[:, index * 2])
        cell = Cell(rawresult[real, index * 2],
                    rawresult[real, index * 2 + 1])
        celllist.append(cell)

# Print all the original data to the report
title = np.empty(setnum * 2, dtype='S20')
title[::2] = 'Concentration:'
title[1::2] = 'Response:'
report.tabletitle(title)
rawresult = np.genfromtxt(filename, delimiter=',')
mask = np.isnan(rawresult)
rawresult = rawresult.astype('string')
rawresult[mask] = ' '
report.table(rawresult)

print 'Select equations.'
print '1. Hill equation'
Equation = check_input('Equation number:', ['1'], 1)
if Equation == 1:
    fitting_curve_hill_equation(filename, celllist, report)
