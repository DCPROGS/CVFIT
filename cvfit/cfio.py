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

