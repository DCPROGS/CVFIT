
from scipy import stats
import numpy as np

class Hill(object):
    def __init__(self, eqname, pars=None):
        """
        pars = [Ymin, Ymax, EC50, nH]
        """
        self.eqname = eqname
        self.pars = pars
        self.fixed = []
        self.names = ['Ymin', 'Ymax', 'EC50', 'nH  ']
        self.data = None
        self.guess = None
        self._theta = None
        
    def equation(self, conc, coeff):
        '''
        The hill equation.
        '''
        return (coeff[0] + (coeff[1] * (conc / coeff[2]) ** coeff[3]) / 
            (1 + (conc / coeff[2]) ** coeff[3]))
            
    def to_fit(self, theta, conc):
        #for each in np.nonzero(self.fixed)[0]:   
        #    theta = np.insert(theta, each, self.pars[each])
        self._set_theta(theta)
        return self.equation(conc, self.pars)
    
    def _set_theta(self, theta):
        for each in np.nonzero(self.fixed)[0]:   
            theta = np.insert(theta, each, self.pars[each])
        self.pars = theta
    def _get_theta(self):
        return self.pars[np.nonzero(np.invert(self.fixed))[0]]
    theta = property(_get_theta, _set_theta)
    
    def propose_guesses(self, set):
        '''
        Calculate the initial guesses for fitting with Hill equation.
        '''
        
        #if self.Component == 1:
        self.guess = np.empty(4)
        slope, intercept, r_value, p_value, std_err = stats.linregress(set.X, set.Y)
        if slope > 0:
            self.trend = 1
        else:
            self.trend = -1
        if self.trend == 1: # Response increases with concentration
            # Determine Y(0)
            if self.fixed[0]:
                self.guess[0] = 0
            else:
                self.guess[0] = np.mean(set.Y[set.X == set.X[0]])
            if self.fixed[1]:
                self.guess[1] = 1
            else:
                # Determine Ymax
                self.guess[1] = np.mean(set.Y[set.X == set.X[-1]]) - self.guess[0]
        else: # Response decreases with concentration
            # Determine Y(0)
            self.guess[0] = np.mean(set.Y[set.X == set.X[-1]])
            # Determine Ymin
            self.guess[1] = np.mean(set.Y[set.X == set.X[0]]) - self.guess[0]
        # Determine Kr
        self.guess[2] = 10 ** ((np.log10(set.X[0]) + np.log10(set.X[-1])) / 2)
        # Determine nH  
        LinRegressX = np.log10(set.X[set.Y < np.amax(set.Y)]) - np.log10(self.guess[2])
        ratio = set.Y[set.Y < np.amax(set.Y)] / np.amax(set.Y)
        LinRegressY = np.log10(ratio / (1 - ratio))
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            LinRegressX, LinRegressY)
        self.guess[3] = slope

#        elif self.Component == 2:
#            print 'Two Components fitting is not completed.'
#            sys.exit(0)
        self.pars = self.guess.copy()
        #return self.guess

