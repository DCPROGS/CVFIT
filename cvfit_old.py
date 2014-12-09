#! /usr/bin/python
"""A general purpose curve fitting software."""
import os
import numpy as np

from hillequation import Cell 
from toolkit import Report, check_input
from hillequation import fitting_curve_hill_equation


print 'Please type in the loaction of the csv file'
filename = raw_input('filename:').strip()

#filename = "Example/Example.csv"

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
