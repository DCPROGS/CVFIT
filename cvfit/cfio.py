import os
import numpy as np
import data

def ask_for_file():
    print 'Please type in the loaction of the csv file'
    filename = raw_input('filename:').strip()
    while (os.path.exists(filename) is False) or (filename[-4:] != '.csv'):
        if os.path.exists(filename) is False:
            print filename, 'not found. Please type in again.'
        elif filename[-4:] != '.csv':
            print filename, 'is not a CSV file. Please type in again.'
        filename = raw_input('filename:').strip()
    return filename

def string_estimates(func, aproxSD):
    j = 0
    str = ''
    for i in range(len(func.names)):
        str += '\nParameter {0:d}: {1}\t= {2:.6g}\t'.format(i+1, func.names[i], func.pars[i])
        if not func.fixed[i]:
            str += 'Approx SD = {0:.6g}'.format(aproxSD[j])
            j += 1
        else:
            str += '(fixed)'
    return str

def string_liklimits(func, limits):
    j = 0
    str = ''
    for i in range(len(func.names)):
        str += '\nParameter {0:d}: {1}\t= {2:.6g}'.format(i+1, func.names[i], func.pars[i])
        if not func.fixed[i]:
            str += '\tLOWER = {0:.6g}'.format(limits[j][0])
            str += '\tUPPER = {0:.6g}'.format(limits[j][1])
            j += 1
        else:
            str += '\t(fixed)'
    return str

def read_sets_from_csv(filename, col=2):
    rawresult = np.genfromtxt(filename, delimiter=',')
    setnum = rawresult.shape[1] / col
    setlist = []
    for i in range(setnum):
        real = np.isfinite(rawresult[:, i * col])
        set = data.XYDataSet()
        if col == 2:
            set.from_columns(rawresult[real, i * col],
               rawresult[real, i * col + 1],
               rawresult[real, i * col] * 0)
        elif col == 3:
            set.from_columns(rawresult[real, i * col],
               rawresult[real, i * col + 1],
               rawresult[real, i * col + 2],)
        set.title = 'Set ' + str(i+1)
        setlist.append(set)
    return setlist