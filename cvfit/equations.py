import math
from scipy import stats
import numpy as np

class Equation(object):
    def __init__(self):
        """        """
        self.eqname = None
        self.ncomp = 1
        self.pars = None
        self.fixed = []
        self.names = []
        self.data = None
        self.guess = None
        self._theta = None
        self.normalised = False
        
    def equation(self, x, coeff):
        ''' '''
        pass
            
    def to_fit(self, theta, x):
        self._set_theta(theta)
        return self.equation(x, self.pars)
    
    def _set_theta(self, theta):
        for each in np.nonzero(self.fixed)[0]:   
            theta = np.insert(theta, each, self.pars[each])
        self.pars = theta
    def _get_theta(self):
        return self.pars[np.nonzero(np.invert(self.fixed))[0]]
    theta = property(_get_theta, _set_theta)
    

class Linear(Equation):
    def __init__(self, eqname, pars=None):
        """
        pars = [a, b]
        """
        self.eqname = eqname
        self.ncomp = 1
        self.pars = pars
        self.fixed = [False, False]
        self.names = ['a', 'b']
        
    def equation(self, x, coeff):
        '''
        The linear equation.
        '''
        return coeff[0] + coeff[1] * x 
   
    def propose_guesses(self, data):
        '''
        Calculate the initial guesses for fitting with Linear equation.
        '''
        #if self.Component == 1:
        slope, intercept, r, p, stderr = stats.linregress(data.X, data.Y)
        self.guess = np.array([intercept, slope])
        self.pars = self.guess.copy()


class Hill(Equation):
    def __init__(self, eqname, pars=None):
        """
        pars = [Ymin, Ymax, EC50, nH]
        """
        self.eqname = eqname
        self.ncomp = 1
        self.pars = pars
        if eqname == 'Hill':
            self.fixed = [True, False, False, False]
            self.names = ['Ymin', 'Ymax', 'EC50', 'nH  ']
        if eqname == 'Langmuir':
            self.fixed = [True, False, False, True]
            self.names = ['Ymin', 'Ymax', 'EC50']
        self.normalised = False
        
    def equation(self, conc, coeff):
        '''
        The hill equation.
        '''
        return (coeff[0] + (coeff[1] * (conc / coeff[2]) ** coeff[3]) / 
            (1 + (conc / coeff[2]) ** coeff[3]))
            
    def normalise(self, data):
        '''
        Nomalise Y to the fitted maximum.
        '''
        # Nomalise the coefficients by fixing the Y(0) and Ymax
        self.normpars = self.pars.copy()
        self.normpars[0], self.normpars[1] = 0, 1
        #data.normY = (data.Y - self.pars[0]) / self.pars[1]
        data.normY = data.Y / self.pars[1]
        data.normS = data.S / self.pars[1]
        self.normalised = True
    
    def propose_guesses(self, data):
        '''
        Calculate the initial guesses for fitting with Hill equation.
        '''
        #if self.Component == 1:
        self.guess = np.empty(4)
        if data.increase: # Response increases with concentration
            # Determine Y(0)
            if self.fixed[0]:
                self.guess[0] = 0
            else:
                self.guess[0] = np.mean(data.Y[data.X == data.X[0]])
            if self.fixed[1]:
                self.guess[1] = 1
            else:
                # Determine Ymax
                self.guess[1] = np.mean(data.Y[data.X == data.X[-1]]) - self.guess[0]
        else: # Response decreases with concentration
            # Determine Y(0)
            if self.fixed[0]:
                self.guess[0] = 0
            else:
                self.guess[0] = np.mean(data.Y[data.X == data.X[-1]])
            # Determine Ymin
            self.guess[1] = np.mean(data.Y[data.X == data.X[0]]) - self.guess[0]
        # Determine Kr
        self.guess[2] = 10 ** ((np.log10(data.X[0]) + np.log10(data.X[-1])) / 2)
        # Determine nH  
        LinRegressX = np.log10(data.X[data.Y < np.amax(data.Y)]) - np.log10(self.guess[2])
        ratio = data.Y[data.Y < np.amax(data.Y)] / np.amax(data.Y)
        LinRegressY = np.log10(ratio / (1 - ratio))
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            LinRegressX, LinRegressY)
        if math.isnan(slope):
            self.guess[3] = 1.0 if data.increase else -1.0
        else:
            self.guess[3] = slope
        if self.eqname == 'Langmuir':
            self.guess[3] = slope / math.fabs(slope)
#        elif self.Component == 2:
#            print 'Two Components fitting is not completed.'
#            sys.exit(0)
        self.pars = self.guess.copy()

