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

def read_sets_from_csv(filename):
    rawresult = np.genfromtxt(filename, delimiter=',')
    setnum = rawresult.shape[1] / 2
    setlist = []
    for i in range(setnum):
        real = np.isfinite(rawresult[:, i * 2])
        set = data.DataSet()
        set.from_columns(rawresult[real, i * 2],
           rawresult[real, i * 2 + 1])
        setlist.append(set)
    return setlist