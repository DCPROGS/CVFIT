
from scipy import stats
import numpy as np

class Hill(object):
    def __init__(self, pars=None):
        """
        pars = [Ymin, Ymax, EC50, nH]
        """
        self.pars = pars
        self.fixed = []
        self.names = ['Ymin', 'Ymax', 'EC50', 'nH  ']
        self.data = None
        self.guess = None
    
    def equation(self, conc):
        '''
        The hill equation.
        '''
        return (self.pars[0] + (self.pars[1] * (conc / self.pars[2]) ** self.pars[3]) / 
            (1 + (conc / self.pars[2]) ** self.pars[3]))
            
    def to_fit(self, theta, conc):
        for each in np.nonzero(self.fixed)[0]:   
            theta = np.insert(theta, each, self.pars[each])
        self.pars = theta
        return self.equation(conc)
    
    def get_theta(self):
        return self.pars[np.nonzero(np.invert(self.fixed))[0]]
    
    def guesses(self, set):
        '''
        Calculate the initial guesses for fitting with Hill equation.
        '''
        
        #if self.Component == 1:
        guess = np.empty(4)
        slope, intercept, r_value, p_value, std_err = stats.linregress(set.X, set.Y)
        if slope > 0:
            self.trend = 1
        else:
            self.trend = -1
        if self.trend == 1: # Response increases with concentration
            # Determine Y(0)
            if self.fixed[0]:
                guess[0] = 0
            else:
                guess[0] = np.mean(set.Y[set.X == set.X[0]])
            if self.fixed[1]:
                guess[1] = 1
            else:
                # Determine Ymax
                guess[1] = np.mean(set.Y[set.X == set.X[-1]]) - guess[0]
        else: # Response decreases with concentration
            # Determine Y(0)
            guess[0] = np.mean(set.Y[set.X == set.X[-1]])
            # Determine Ymin
            guess[1] = np.mean(set.Y[set.X == set.X[0]]) - guess[0]
        # Determine Kr
        guess[2] = 10 ** ((np.log10(set.X[0]) + np.log10(set.X[-1])) / 2)
        # Determine nH  
        LinRegressX = np.log10(set.X[set.Y < np.amax(set.Y)]) - np.log10(guess[2])
        ratio = set.Y[set.Y < np.amax(set.Y)] / np.amax(set.Y)
        LinRegressY = np.log10(ratio / (1 - ratio))
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            LinRegressX, LinRegressY)
        guess[3] = slope

#        elif self.Component == 2:
#            print 'Two Components fitting is not completed.'
#            sys.exit(0)
        self. pars = guess.copy()
        return guess

