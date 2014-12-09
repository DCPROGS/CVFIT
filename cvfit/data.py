import numpy as np

class XYDataSet(object):
    """
    """
    def __init__(self):
        self.X = []
        self.Y = []
        self.SE = []
        self.W = []
        self.title = None
        self._weightmode = None

    def from_columns(self, X, Y, SE=None):
        # Keep track of the original data
        self.Xoriginal, self.Yoriginal, self.SEoriginal = X, Y, SE
        # Sort the data according to concentration
        temp = np.vstack((self.Xoriginal, self.Yoriginal, self.SEoriginal))
        temp.sort()
        self.X = temp[0]
        self.Y = temp[1]
        self.SE = temp[2]
        self.W = np.ones((len(self.X)))

    def size(self):
        return len(self.X)
    
    def _set_weightmode(self, weightmode):
        if weightmode == 1:
            pass
        elif weightmode == 2:
            self.W = np.W / np.power(np.array(self.SE), 2)
        self._weightmode = weightmode 
    def _get_weightmode(self):
        return self._weightmode
    weightmode = property(_get_weightmode, _set_weightmode)

    def __str__(self):
        str = "\n"
        for i in range(len(self.X)):
            str += "{0:.6g}\t{1:.6g}\n".format(self.X[i], self.Y[i], self.SE[i])
        return str