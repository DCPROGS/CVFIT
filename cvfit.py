#! /usr/bin/python

__author__="remislp"
__date__ ="$29-Nov-2014 15:40:02$"

from cvfit.hill import Hill
from cvfit import fitting

if __name__ == "__main__":
    
    sets = fitting.load_data(example=True)
    settings = fitting.general_settings()
    fixed = [True, False, False, False]
    for set in sets:
        hill = Hill()
        hill.fixed = fixed
        fitting.session(set, hill)
