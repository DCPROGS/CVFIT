__author__="remis"
__date__ ="$23-Feb-2010 10:29:31$"

import sys
from PySide.QtGui import *
from PySide.QtCore import *

from cvfit.fitting import SingleFitSession

class EquationDlg(QDialog):
    def __init__(self, data, eqtype='Hill', eqname='Hill', output=sys.stdout, parent=None):
        super(EquationDlg, self).__init__(parent)
        self.setWindowTitle(eqname + ' equation settings...')
        if eqtype == 'Hill':
            from cvfit.hill import Hill as EQ
        self.equations = []
        self.data = data
        self.log = output
        
        layoutMain = QVBoxLayout()
        self.layouts = []
        for i in range(len(self.data)):
            equation = EQ(eqname)
            layoutMain.addWidget(QLabel(self.data[i].title))
            layout = QHBoxLayout()
            layout.addWidget(QLabel("Number of components:"))
            compSpinBox = QSpinBox()
            compSpinBox.setAlignment(Qt.AlignRight|Qt.AlignVCenter)
            compSpinBox.setRange(1, 2)
            compSpinBox.setValue(equation.ncomp)
            self.connect(compSpinBox, SIGNAL("valueChanged(int)"), self.on_setting_changed)
            layout.addWidget(compSpinBox)
            y0fixedChB = QCheckBox("&Fix Y(0) or Y(inf) = 0?")
            y0fixedChB.setChecked(True)
            self.connect(y0fixedChB, SIGNAL("stateChanged(int)"), self.on_setting_changed)
            layout.addWidget(y0fixedChB)
            normChB = QCheckBox("&Data normalised?")
            normChB.setChecked(equation.normalised)
            self.connect(normChB, SIGNAL("stateChanged(int)"), self.on_setting_changed)
            layout.addWidget(normChB)
            self.equations.append(equation)
            self.layouts.append(layout)
            layoutMain.addLayout(layout)

        buttonBox = QDialogButtonBox(QDialogButtonBox.Ok|QDialogButtonBox.Cancel)
        buttonBox.button(QDialogButtonBox.Ok).setDefault(True)
        self.connect(buttonBox, SIGNAL("accepted()"),
             self, SLOT("accept()"))
        self.connect(buttonBox, SIGNAL("rejected()"),
             self, SLOT("reject()"))
        buttonLayout = QHBoxLayout()
        buttonLayout.addStretch()
        buttonLayout.addWidget(buttonBox)
        layoutMain.addLayout(buttonLayout)
        self.setLayout(layoutMain)

    def on_setting_changed(self):
        for i in range(len(self.data)):
            self.equations[i].ncomp = self.layouts[i].itemAt(0).widget().value()
            if self.layouts[i].itemAt(1).widget().isChecked():
                self.equations[i].fixed[0] = self.layouts[i].itemAt(1).widget().isChecked()
                self.equations[i].pars[0] = 0.0
            self.equations[i].normalised = self.layouts[i].itemAt(2).widget().isChecked()
        
    def return_equation(self):
        fsessions = []
        for i in range(len(self.data)):
            fsessions.append(SingleFitSession(self.data[i], self.equations[i], self.log))
        return fsessions

class GuessDlg(QDialog):
    def __init__(self, fs, parent=None):
        super(GuessDlg, self).__init__(parent)
        self.fs = fs
        layoutMain = QVBoxLayout()

        self.layouts = []
        for i in range(len(self.fs)):
            layoutMain.addWidget(QLabel(fs[i].data.title))
            layout = QHBoxLayout()
            for j in range(len(self.fs[i].eq.names)):
                layout.addWidget(QLabel(self.fs[i].eq.names[j]))
                pared = QLineEdit(unicode(self.fs[i].eq.guess[j]))
                pared.setMaxLength(6)
                self.connect(pared, SIGNAL("editingFinished()"), self.on_guess_changed)
                layout.addWidget(pared)
                parfix = QCheckBox("Fixed?")
                parfix.setChecked(self.fs[i].eq.fixed[j])
                self.connect(parfix, SIGNAL("stateChanged(int)"), self.on_guess_changed)
                layout.addWidget(parfix)
            self.layouts.append(layout)
            layoutMain.addLayout(layout)

        buttonBox = QDialogButtonBox(QDialogButtonBox.Ok|QDialogButtonBox.Cancel)
        buttonBox.button(QDialogButtonBox.Ok).setDefault(True)
        self.connect(buttonBox, SIGNAL("accepted()"),
             self, SLOT("accept()"))
        self.connect(buttonBox, SIGNAL("rejected()"),
             self, SLOT("reject()"))
        layoutMain.addWidget(buttonBox)
        self.setLayout(layoutMain)
        self.setWindowTitle('Edit guesses...')

    def on_guess_changed(self):
        for i in range(len(self.fs)):
            for j in range(len(self.fs[i].eq.names)):
                self.fs[i].eq.guess[j] = float(self.layouts[i].itemAt(j*3+1).widget().text())
                self.fs[i].eq.fixed[j] = self.layouts[i].itemAt(j*3+2).widget().isChecked()

    def return_guesses(self):
        return self.fs
