__author__="remis"
__date__ ="$23-Feb-2010 10:29:31$"

import sys
import csv
import xlrd
from PyQt4.QtGui import *
from PyQt4.QtCore import *

from cvfit import data
from cvfit.fitting import SingleFitSession, MultipleFitSession

class WarningDlg(QDialog):
    def __init__(self, message, parent=None):
        super(WarningDlg, self).__init__(parent)
        self.setWindowTitle('Warning...')
 
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel(message))
        
        buttonBox = QDialogButtonBox(QDialogButtonBox.Ok|QDialogButtonBox.Cancel)
        buttonBox.button(QDialogButtonBox.Ok).setDefault(True)
        self.connect(buttonBox, SIGNAL("accepted()"),
             self, SLOT("accept()"))
        self.connect(buttonBox, SIGNAL("rejected()"),
             self, SLOT("reject()"))
        layout.addWidget(buttonBox)


class LoadDataDlg(QDialog):
    def __init__(self, filename, sheet=None, parent=None):
        super(LoadDataDlg, self).__init__(parent)
        self.setWindowTitle('Check data format...')
        self.resize(600, 300)
        layout = QVBoxLayout(self)
        
        self.filename = filename
        self.sheet = sheet
        self.type = filename.split('.')[-1]
        self.col = 2
        self.row = 0
        self.weight = 1
        
        textBox = QTextBrowser()
        layout.addWidget(textBox)
        
        if self.type == 'xls' or self.type == 'xlsx':
            wb = xlrd.open_workbook(filename)
            s = wb.sheet_by_index(sheet)
            txtdata = ''
            for curr_row in range(s.nrows):
                 row = s.row(curr_row)
                 rowstr = ''
                 for cell in row:
                     rowstr += (str(cell.value) + '\t')
                 txtdata += (rowstr + '\n')
            textBox.append(txtdata)
        else:
            f = open(filename, 'rU')
            txtdata = list(csv.reader(f, dialect=csv.excel_tab))
            f.close()
            for row in txtdata:
                textBox.append(''.join(row))
        
        layout1 = QHBoxLayout()
        layout1.addWidget(QLabel('How many header lines to skip?'))
        self.rowSB = QSpinBox()
        self.rowSB.setValue(0)
        self.rowSB.setRange(0, 100)
        self.rowSB.valueChanged.connect(self.on_changed)
        layout1.addWidget(self.rowSB)
        
        layout2 = QHBoxLayout()
        layout2.addWidget(QLabel('How many columns in each set (2 or 3)?'))
        self.colSB = QSpinBox()
        self.colSB.setValue(2)
        self.colSB.setRange(2, 3)
        self.colSB.valueChanged.connect(self.on_changed)
        layout2.addWidget(self.colSB)
        layout.addLayout(layout1)
        layout.addLayout(layout2)
        
        layout.addWidget(QLabel("Please select the weighting method now:"))
        self.group1 = QButtonGroup(self)
        self.rb1 = QRadioButton("Weights constant; errors from residuals (Default).")
        self.rb1.setChecked(True)
        self.group1.addButton(self.rb1)
        layout.addWidget(self.rb1)
        self.rb2 = QRadioButton("Weights from specified s(Y); errors from weights.")
        self.group1.addButton(self.rb2)
        layout.addWidget(self.rb2)
        self.group1.buttonClicked.connect(self.on_changed)
        
        buttonBox = QDialogButtonBox(QDialogButtonBox.Ok|QDialogButtonBox.Cancel)
        buttonBox.button(QDialogButtonBox.Ok).setDefault(True)
        self.connect(buttonBox, SIGNAL("accepted()"),
             self, SLOT("accept()"))
        self.connect(buttonBox, SIGNAL("rejected()"),
             self, SLOT("reject()"))
        layout.addWidget(buttonBox)     
             
    def on_changed(self):
        if self.group1.checkedId() == -2:
            self.weight = 1
        elif self.group1.checkedId() == -3:
            self.weight = 2
        self.col = self.colSB.value()
        self.row = self.rowSB.value()
        
    def return_data(self):
        if self.type == 'xls' or self.type == 'xlsx':
            return data.read_sets_from_Excel(self.filename, self.col, self.row, 
                self.sheet)
        else:
            return data.read_sets_from_csv(self.filename, self.type, col=self.col,
                header=self.row, namesin=False, weight=self.weight)
            
class ExcelSheetDlg(QDialog):
    """
    Dialog to choose Excel sheet to load.
    """
    def __init__(self, sheetlist, parent=None):
        super(ExcelSheetDlg, self).__init__(parent)

        self.sheet = ''
        self.List = QListWidget()
        self.List.addItems(sheetlist)
    
        self.connect(self.List,
            SIGNAL("itemSelectionChanged()"),
            self.sheetSelected)

        buttonBox = QDialogButtonBox(QDialogButtonBox.Ok|
            QDialogButtonBox.Cancel)
        self.connect(buttonBox, SIGNAL("accepted()"),
            self, SLOT("accept()"))
        self.connect(buttonBox, SIGNAL("rejected()"),
            self, SLOT("reject()"))

        layout1 = QHBoxLayout()
        layout1.addWidget(self.List)
        layout2 = QVBoxLayout()
        layout2.addLayout(layout1)
        layout2.addWidget(buttonBox)

        self.setLayout(layout2)
        self.resize(200, 300)
        self.setWindowTitle("Choose Excel sheet to load...")

    def sheetSelected(self):
        """
        Get selected sheet name.
        """
        self.sheet = self.List.currentRow()

    def returnSheet(self):
        """
        Return selected sheet name.
        """
        return self.sheet
        

class EquationDlg(QDialog):
    def __init__(self, data, eqtype='Hill', eqname='Hill', output=sys.stdout, parent=None):
        super(EquationDlg, self).__init__(parent)
        self.setWindowTitle(eqname + ' equation settings...')
        if eqtype == 'Hill':
            from cvfit.equations import Hill as EQ
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
        fits = MultipleFitSession(self.log)
        for i in range(len(self.data)):
            fits.add(SingleFitSession(self.data[i], self.equations[i], self.log))
        return fits

class GuessDlg(QDialog):
    def __init__(self, fs, parent=None):
        super(GuessDlg, self).__init__(parent)
        self.fs = fs
        layoutMain = QVBoxLayout()

        self.layouts = []
        for i in range(len(self.fs.list)):
            layoutMain.addWidget(QLabel(fs.list[i].data.title))
            layout = QHBoxLayout()
            for j in range(len(self.fs.list[i].eq.names)):
                layout.addWidget(QLabel(self.fs.list[i].eq.names[j]))
                pared = QLineEdit(str(self.fs.list[i].eq.guess[j]))
                pared.setMaxLength(6)
                self.connect(pared, SIGNAL("editingFinished()"), self.on_guess_changed)
                layout.addWidget(pared)
                parfix = QCheckBox("Fixed?")
                parfix.setChecked(self.fs.list[i].eq.fixed[j])
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
        for i in range(len(self.fs.list)):
            for j in range(len(self.fs.list[i].eq.names)):
                self.fs.list[i].eq.guess[j] = float(self.layouts[i].itemAt(j*3+1).widget().text())
                self.fs.list[i].eq.fixed[j] = self.layouts[i].itemAt(j*3+2).widget().isChecked()

    def return_guesses(self):
        for i in range(len(self.fs.list)):
            self.fs.list[i].eq.pars = self.fs.list[i].eq.guess
        return self.fs
