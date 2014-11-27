"""
Equations available to fit data in CVFIT.
"""
import copy
import numpy as np
def SSD(gars, func, data, notfixed, fixval):
    """
    Calculate sum of squared deviations.
    """
   # print '\n data \n', data 
    S = 0.0
   # if type(data) is np.ndarray:
   #     S = 2
   # else:
    for point in data.points:
        S += (point.y - func(gars, point.x, notfixed, fixval)) ** 2
    return S
    
def hill_equation(gars, conc, notfixed, fixval):
    '''
    The hill equation.
    '''
    #print '\n gars, notfixed \n', gars, notfixed
    for i in range(0,len(gars)):
        if not notfixed[i]:
            gars[i] = copy.deepcopy(fixval[i])
    ymin, ymax, ec50, nH = gars
    return ymin + (ymax * (conc / ec50) ** nH) / (1 + (conc / ec50) ** nH)
    



