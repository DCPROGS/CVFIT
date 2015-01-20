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

def check_input(text, accept, default):
    '''
    Check if the input is in acceptable range or not
    If not, ask to key in another value
    '''

    inputnumber = raw_input(text)
    if inputnumber:
        while inputnumber not in accept:
            print text
            inputnumber = raw_input(text)
    else:
        inputnumber = default

    return int(inputnumber)

def string_estimates(func, aproxSD, CVs):
    j = 0
    str = ''
    for i in range(len(func.names)):
        str += '\nParameter {0:d}: {1}\t= {2:.6g}\t'.format(i+1, func.names[i], func.pars[i])
        if not func.fixed[i]:
            str += 'Approx SD = {0:.6g}\t'.format(aproxSD[j])
            str += 'CV = {0:.1f}'.format(CVs[j])
            j += 1
        else:
            str += '(fixed)'
    if np.any(CVs > 33):
        str += "\nWARNING: SOME PARAMETERS POORLY DEFINED (CV > 33%); try different guesses"
    return str

def string_liklimits(func, limits):
    j = 0
    str = ''
    for i in range(len(func.names)):
        str += '\nParameter {0:d}: {1}\t= {2:.6g}'.format(i+1, func.names[i], func.pars[i])
        if not func.fixed[i]:
            try:
                str += '\tLOWER = {0:.6g}'.format(limits[j][0])
            except:
                str += '\tLOWER limit not found'
            try:
                str += '\tUPPER = {0:.6g}'.format(limits[j][1])
            except:
                str += '\tUPPER limit not found'
            j += 1
        else:
            str += '\t(fixed)'
    return str

def read_sets_from_csv(filename, type, col=2, header=0, namesin=False, weight=1):
    """
    Read data from a txt or csv file in which X and Y (and SD) are in columns.
    """
    if type == 'csv':
        delimit = ','
    elif type == 'txt':
        delimit = '\t'
    if namesin:
        f = open(filename, 'r')
        nameline = f.readline().strip("\n")
        names = nameline.split(delimit)
        f.close()

    rawresult = np.genfromtxt(filename, delimiter=delimit,skip_header=header)
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
        set.weightmode = weight
        setlist.append(set)
    return setlist

