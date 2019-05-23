import xlrd
import numpy as np
import pandas as pd
from scipy import stats

class XYDataSet(object):
    """
    """
    def __init__(self):
        self.X = []
        self.Y = []
        self.S = [] # Standard deviation of the mean
        self.W = [] # Weight
        self.title = '*no*name*'
        self._weightmode = 1
        self.increase = None

    def from_columns(self, X, Y, S=None):
        self.X, self.Y = np.array(X), np.array(Y)
        if S is None:
            self.S = np.zeros((len(self.X)))
        else:
            self.S = np.array(S)
        if self.S.any() == 0:
            self.W = np.ones((len(self.X)))
        else:
            self.W = 1.0 / np.power(np.array(self.S), 2)
            self._weightmode = 2
        self.sort()
        self._set_trend()
        
    def pool(self, X, Y, S):
        self.X, self.Y = np.hstack((self.X, X)), np.hstack((self.Y, Y))
        self.S, self.W = np.hstack((self.S, S)), np.ones((len(self.X)))
        self.sort()
        self._set_trend()
        self.title = 'Pooled data'
        
    def sort(self):
        # Sort the data according to concentration
        
        idx = np.argsort(self.X)
        self.X, self.Y = self.X[idx], self.Y[idx]
        self.S, self.W = self.S[idx], self.W[idx]
        
    def average_pooled(self):
        '''
        Calculate pooled data average response for each concentration and stardard
        deviation of the mean.
        '''

        self.avX = np.unique(self.X)
        avY, avS = [], []
        for con in self.avX:
            Y = self.Y[self.X == con]
            if len(Y) > 1:
                avY.append(np.mean(Y))
                avS.append(np.std(Y) / np.sqrt(len(Y)))
            else:
                avY.append(Y[0])
                avS.append(0.0)
        self.avY = np.array(avY)
        self.avS = np.array(avS)

    def size(self):
        return len(self.X)
    def Xmin(self):
        return np.min(self.X)
    def Xmax(self):
        return np.max(self.X)
    def Ymin(self):
        return np.min(self.Y)
    def Ymax(self):
        return np.max(self.Y)
    
    def _set_trend(self):
        slope, intercept, r_value, p_value, std_err = stats.linregress(self.X, self.Y)
        if slope > 0:
            self.increase = True
        else:
            self.increase = False
    
    def _set_weightmode(self, weightmode):
        """
        Set weights of observations. 
        Weightmode 1:  Weights equal for all points.  In this case the weights
        are set to 1.0, and give no information about the precision of the
        measurements.  For the calculation of errors in the parameter estimates,
        the scatter is estimated from the extent to which the data points fail
        to lie exactly on the fitted curve, i.e. from the 'residuals'  (so this
        method assumes that an equation has been chosen that really describes
        the data).
        Weightmode 2:  Weights are calculated as inverse of the square of the
        standard deviation (the variance) supplied during data loading.  This is
        the best method when a legitimate estimate of the precision of each
        point is known.  Usually the estimates will come from replicate
        observations at each X value, which are averaged, and the standard
        deviation of that mean used to calculate the weights.
        Weightmodes 3 and 4 not implemented yet.
        Weightmode 3:  Weights are calculated as in Weightmode 2; the standard
        deviation (the variance) is calculated from Y repeats at the same X.
        Weightmode 4:  Assume that all individual observations have equal
        precision and specify n, the number of values in the average, as the
        relative weight.  Only relative values are given so they need contain no
        information about the size of the experimental error (estimate errors
        are calculated from residuals as for method 1).  One case in which this
        may be useful is when the data points are means, with different numbers
        of observations averaged for each.  
        """
        
        self._weightmode = weightmode 
        if weightmode == 1: 
            self.W = np.ones((len(self.X)))
        elif weightmode == 2:
            if self.S.any() == 0:
                self.W = np.ones((len(self.X)))
                self._weightmode = 1
                print('data: WARNING : some SD equal to 0;' + 
                    'cannot be used for weights; reverting to weightmode = 1')
            else:
                self.W = 1.0 / np.power(np.array(self.S), 2)
    def _get_weightmode(self):
        return self._weightmode
    weightmode = property(_get_weightmode, _set_weightmode)

    def __str__(self):
        str = "\nX\tY\ts(Y)\tweight\n"
        for i in range(len(self.X)):
            str += "{0:.6g}\t{1:.6g}\t{2:.6g}\t{3:.6g}\n".format(self.X[i], self.Y[i], self.S[i], self.W[i])
        return str

def read_single_columns_from_Excel(fname, sheet):
    """
    Read Excel file sheet where each column contains multiple measurements for one condition
    (eg concentration).
    """
    xl = pd.ExcelFile(fname)
    df = xl.parse(sheet)
    names = df.columns.tolist()
    setlist = []
    X, Y, S = [], [], []
    for i in range(len(names)):
        temp = names[i].split()
        x, unit = float(temp[0]), temp[1]
        y = df.iloc[:, i].dropna().values.tolist()
        Y += y
        X += [x] * len(y)
        S += [0] * len(y)
    set = XYDataSet()
    set.from_columns(np.array(X), np.array(Y), np.array(S))
    set.title = 'Set ' + str(sheet)
    set.weightmode = 1
    setlist.append(set)
    return setlist
    
def read_sets_from_Excel(fname, set_col=0, line_skip=0, sheet=0):
    """
    Read Excel file sheet to load data.
    """
    
    wb = xlrd.open_workbook(fname)
    s = wb.sheet_by_index(sheet)
    setlist = []
    for i in range(int(s.ncols / set_col)):
        X, Y, S = [], [], []
        if set_col == 2:
            for cell1, cell2 in zip(s.col_slice(i*set_col, start_rowx=line_skip), s.col_slice(i*set_col+1, start_rowx=line_skip)):
                if cell1.ctype not in (xlrd.XL_CELL_EMPTY, xlrd.XL_CELL_BLANK):
                    X.append(cell1.value), Y.append(cell2.value), S.append(0)
        if set_col == 3:
            for cell1, cell2, cell3 in zip(s.col_slice(i*set_col, start_rowx=line_skip),
                s.col_slice(i*set_col+1, start_rowx=line_skip),
                s.col_slice(i*set_col+2, start_rowx=line_skip)):
                if cell1.ctype not in (xlrd.XL_CELL_EMPTY, xlrd.XL_CELL_BLANK):
                    X.append(cell1.value), Y.append(cell2.value), S.append(cell3.value)
        set = XYDataSet()
        set.from_columns(np.array(X), np.array(Y), np.array(S))
        set.title = 'Set ' + str(i+1)
        set.weightmode = 1
        setlist.append(set)
    return setlist


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
        set = XYDataSet()
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


def ask_for_file():
    print('Please type in the loaction of the csv file')
    filename = raw_input('filename:').strip()
    while (os.path.exists(filename) is False) or (filename[-4:] != '.csv'):
        if os.path.exists(filename) is False:
            print(filename, 'not found. Please type in again.')
        elif filename[-4:] != '.csv':
            print (filename, 'is not a CSV file. Please type in again.')
        filename = raw_input('filename:').strip()
    return filename