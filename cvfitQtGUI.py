#! /usr/bin/python

"""
Load set(s) of data from an Excel, csv or text file.
From a checkable list select set(s) of data to fit.
Plot the data.
Choose a function to fit the data points.
Fit the function to the data using ...
Calculate approximate SD and likelihood intervals.
Display the fit on the data plot.
Generate HTML report.

Supported data input file formats:
1. Excel .xls and .xlsx files.
2. CSV file.
3. TXT file.

"""

__author__="remislp"
__date__ ="$07-Feb-2010 16:23:20$"

import sys
try:
    from PySide.QtGui import *
    from PySide.QtCore import *
except:
    raise ImportError("pyqt module is missing")

import cvfit.myqtlib.blocks as blocks

class CVFITgui(QMainWindow):
    def __init__(self, parent=None):
        super(CVFITgui, self).__init__(parent)
        self.resize(1000, 700)     # wide, high in px
        self.mainFrame = QWidget()
        self.setWindowTitle("CVFIT: Curve fitting")
        
        self.log = None
        self.data = None
        self.fits = None
        self.pooledfit = None
        self.eqfit = None
        self.eqname = None
        self.eqtype = None
        self.fname = ''
        self.report = None
        
        upperHBox = QHBoxLayout()
        self.datablk = blocks.DataBlock(self)
        upperHBox.addWidget(self.datablk)
        self.eqblk = blocks.EquationBlock(self)
        upperHBox.addWidget(self.eqblk)
        
        lowerHBox = QHBoxLayout()
        self.textblk = blocks.TextBlock(self)
        lowerHBox.addWidget(self.textblk)
        self.plotblk = blocks.PlotBlock(self)
        lowerHBox.addWidget(self.plotblk)
        
        VBox = QVBoxLayout()
        VBox.addLayout(upperHBox)
        VBox.addLayout(lowerHBox)
        self.mainFrame.setLayout(VBox)
        self.setCentralWidget(self.mainFrame)
        
if __name__ == "__main__":
    app = QApplication(sys.argv)
    form = CVFITgui()
    form.show()
    app.exec_()