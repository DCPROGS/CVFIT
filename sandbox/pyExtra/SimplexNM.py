
__author__="Lab"
__date__ ="$02-Sep-2009 13:51:30$"


class Simplex(object):
    'Multidimensional minimisation of the function f(x) where x is a vector \
    in n dimensions, by the downhill simplex method of Nelder and Mead. \
    This version is addopted from Numerical Recipes in Fortran77. \
    Inputs: n- number of free parameters; simp- a matrix [n+1, n], its rows are \
    the vertices of the starting simplex; y- a vector [n+1], whose components \
    are pre-initialised to the values of f(x) evaluated at the n+1 vertices of \
    simp; ftol- the fractional convergence tolerance to be achieved in the f(x) \
    value. Outputs: simp and y, with their values reset to n+1 new points all \
    within ftol of a minimum function value; iter- the number of function \
    evaluations.'



    def __init__(self):

        iter = 0

if __name__ == "__main__":
    print "Hello";