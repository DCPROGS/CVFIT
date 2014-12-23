import numpy as np

class XYDataSet(object):
    """
    """
    def __init__(self):
        self.X = []
        self.Y = []
        self.S = [] # Standard deviation of the mean
        self.W = [] # Weight
        self.title = None
        self._weightmode = 1

    def from_columns(self, X, Y, S=None):
        # Keep track of the original data
        self.Xoriginal, self.Yoriginal, self.Soriginal = X, Y, S
        # Sort the data according to concentration
        temp = np.vstack((self.Xoriginal, self.Yoriginal, self.Soriginal))
        temp.sort()
        self.X = temp[0]
        self.Y = temp[1]
        self.S = temp[2]
        self.W = np.ones((len(self.X)))

    def size(self):
        return len(self.X)
    
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
                self.W = np.array(self.W) / np.power(np.array(self.S), 2)
    def _get_weightmode(self):
        return self._weightmode
    weightmode = property(_get_weightmode, _set_weightmode)

    def __str__(self):
        str = "\nX\tY\ts(Y)\tweight\n"
        for i in range(len(self.X)):
            str += "{0:.6g}\t{1:.6g}\t{2:.6g}\t{3:.6g}\n".format(self.X[i], self.Y[i], self.S[i], self.W[i])
        return str