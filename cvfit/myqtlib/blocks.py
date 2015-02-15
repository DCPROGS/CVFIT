import os
import xlrd
import numpy as np

from PySide.QtGui import *
from PySide.QtCore import *

from cvfit import data
from cvfit import plots
from cvfit.fitting import SingleFitSession
import cvfit.myqtlib.myqtcommon as mqt
import cvfit.myqtlib.dialogs as dialogs
from cvfit.report import Report


class DataBlock(QWidget):
    def __init__(self, parent=None):
        super(DataBlock, self).__init__(parent)
        self.parent = parent
   
        loadButton = QPushButton("&Load new data")
        self.connect(loadButton, SIGNAL('clicked()'), self.on_load)
        dataListView = QListView()
        self.dataListModel = QStandardItemModel()
        dataListView.setModel(self.dataListModel)
        plotButton = QPushButton("&Plot data")
        self.connect(plotButton, SIGNAL('clicked()'), self.on_plot)
        layout = QVBoxLayout(self)
        layout.addWidget(loadButton)
        layout.addWidget(dataListView)
        layout.addWidget(plotButton)
        layout.addStretch(1)
        
    def on_plot(self):
        
        self.parent.data = []
        for row in range(self.dataListModel.rowCount()):
            i = self.dataListModel.index(row, 0)
            checked = (self.dataListModel.data(i, Qt.CheckStateRole) ==
                Qt.Checked)
            if checked:
                self.parent.data.append(self.allsets[row])
        #self.parent.plotblk.on_show()
        plots.plot(self.parent)

    def on_load(self):
        filename, path = QFileDialog.getOpenFileName(self.parent,
                "Open a data file...", ".",
                "Excel files (*.xls *.xlsx);;CSV files (*.csv);" +
                ";TXT files (*.txt);;All files (*.*)")
        self.parent.log.write('\nLoading file: ' + filename)
        type = filename.split('.')[-1]
        
        xlssheet = None
        if type == 'xls' or type == 'xlsx':
            book = xlrd.open_workbook(filename)
            sheets = book.sheet_names()
            dialog = dialogs.ExcelSheetDlg(sheets, self)
            if dialog.exec_():
                xlssheet = dialog.returnSheet()

        dialog = dialogs.LoadDataDlg(filename, sheet=xlssheet)
        if dialog.exec_():
            self.allsets = dialog.return_data()

        #self.allsets = cfio.read_sets_from_csv(filename, col=2)
        self.parent.log.write("Loaded: " + 
            os.path.split(str(filename))[1])
        self.dataListModel.clear()
        for set in self.allsets:
            item = QStandardItem(set.title)
            item.setCheckState(Qt.Checked)
            item.setCheckable(True)
            self.dataListModel.appendRow(item)
            self.parent.log.write('\n'+set.title)
            self.parent.log.write(str(set))

        self.parent.data = self.allsets
        self.parent.fname = filename


class EquationBlock(QWidget):
    def __init__(self, parent=None):
        super(EquationBlock, self).__init__(parent)
        self.parent = parent
        
        eqLabel = QLabel("Please choose equation:                                     ")
        self.eqTitles = [
        '(1) Hill equation(s) (inc. or dec./common K or max)',
        '(2) Langmuir hyperbola(s) plus straight line',
        #'(3) Polynomial (inc. straight line)',
        #'(4) Langmuir hyperbola(s) (inc. or dec.)',
        #'(5) Hill equation(s) plus straight line',
        #'(6) Power function y=ybar*(x/x0)^n (linear on log-log plot)',
        #'(7) Binding inhibition curve (parameter=KB)'
        ''
        ]
        self.eqList = QListWidget()
        for i in range(len(self.eqTitles)):
            item = QListWidgetItem(self.eqTitles[i])
            self.eqList.addItem(item)
            if i == 0:
                item.setSelected(True)

        eqButton = QPushButton("&Load equation")
        self.connect(eqButton, SIGNAL('clicked()'), self.on_equation)
        guessButton = QPushButton("Initial guesses")
        self.connect(guessButton, SIGNAL('clicked()'), self.on_guess)
        fitButton = QPushButton("Fit")
        self.connect(fitButton, SIGNAL('clicked()'), self.on_fit)
        normButton = QPushButton("Normalise and replot")
        self.connect(normButton, SIGNAL('clicked()'), self.on_normalise)
        poolButton = QPushButton("Pool normalised and fit")
        self.connect(poolButton, SIGNAL('clicked()'), self.on_pool)
        
        layout = QVBoxLayout(self)
        layout.addWidget(eqLabel)
        layout.addWidget(self.eqList)
        layout.addWidget(eqButton)
        layout.addWidget(guessButton)
        layout.addWidget(fitButton)
        layout.addWidget(normButton)
        layout.addWidget(poolButton)
        layout.addStretch(1)
        
    def on_pool(self):
        pooldata = data.XYDataSet()
        for session in self.parent.fits:
            pooldata.pool(session.data.X, session.data.normY, session.data.S)
        pooldata.weightmode = 1
        fsession = SingleFitSession(pooldata, self.parent.eqfit, self.parent.log)
        fsession.fit()
        fsession.calculate_errors()
        fsession.data.average_pooled()
        self.parent.pooledfit = fsession

        self.parent.log.write('\n***********************\n**************************')
        self.parent.log.write('\tNormalised and pooled data fit finished')
        self.parent.log.write(fsession.string_estimates())
        self.parent.log.write(fsession.string_liklimits())

        plot_filename = 'fitted_normalised_pooled.png'
        plots.plot(self.parent, plotPooled=True, save_fig_name=plot_filename)
        self.parent.report.title('Pooled data fit finished', 1)
        self.parent.report.paragraph(fsession.string_estimates())
        self.parent.report.paragraph(fsession.string_liklimits())
        self.parent.report.image(plot_filename)
        
    def on_normalise(self):
        for session in self.parent.fits:
            session.eq.normalise(session.data)

        plot_filename = 'all_fitted_normalised_curves.png'
        plots.plot(self.parent, plotNorm=True, save_fig_name=plot_filename)
        self.parent.report.title('Data normalised to the fitted maxima', 1)
        self.parent.report.image(plot_filename)
        
    def on_fit(self):

        self.parent.report = Report(self.parent.fname)
        self.parent.report.title('Original data:', 1)
        self.parent.report.paragraph('Number of datasets loaded: ' + str(len(self.parent.fits)))
        for fit in self.parent.fits:
            self.parent.report.dataset(fit.data.title, str(fit.data))

        progressDlg = QProgressDialog('Fitting data set {0:d}'.format(1),
                                 "Cancel", 1, len(self.parent.fits))
        progressDlg.setWindowTitle('Fitting...')
        i = 1
        for fs in self.parent.fits:
            fs.fit()
            fs.calculate_errors()
            self.parent.log.write('\n*************************************************')
            self.parent.log.write('\t' + fs.data.title + ' fit finished')
            self.parent.log.write(fs.string_estimates())
            self.parent.log.write(fs.string_liklimits())
            self.parent.report.title('{0} fit finished'.format(fs.data.title), 1)
            self.parent.report.paragraph(fs.string_estimates())
            self.parent.report.paragraph(fs.string_liklimits())
            
            progressDlg.setValue(i)
            i += 1
            progressDlg.setLabelText('Fitting data set {0:d}'.format(i))
            if progressDlg.wasCanceled():
                break

        plot_filename = 'all_fittedcurves.png'
        plots.plot(self.parent, plotFit=True, save_fig_name=plot_filename)
        self.parent.report.image(plot_filename)
        
    def on_guess(self):
        plots.plot(self.parent, plotGuesses=True)
        dialog = dialogs.GuessDlg(self.parent.fits)
        if dialog.exec_():
            fs = dialog.return_guesses()
        self.parent.fits = fs
        plots.plot(self.parent, plotGuesses=True)
        
    def on_equation(self):
        row = self.eqList.currentRow()
        if row == 0:
            eqname = 'Hill'
            eqtype = 'Hill'
            from cvfit.hill import Hill as eqfit
        elif row == 1:
            eqname = 'Langmuir'
            eqtype = 'Hill'
            from cvfit.hill import Hill as eqfit
        else:
            self.parent.log.write("This eqation is not implemented yet.")
            self.parent.log.write("Please, choose other equation.")
            
        self.parent.eqname = eqname
        self.parent.eqtype = eqtype
        self.parent.eqfit = eqfit(eqname)
        self.parent.fits = []
        dialog = dialogs.EquationDlg(self.parent.data, eqtype, eqname, 
            self.parent.log)
        if dialog.exec_():
            self.parent.fits = dialog.return_equation()
        
        
class TextBlock(QWidget):
    def __init__(self, parent=None):
        super(TextBlock, self).__init__(parent)
        self.parent = parent
        
        # Prepare text box for printout and set where printout goes
        self.textBox = QTextBrowser()
        self.parent.log = mqt.PrintLog(self.textBox) #, sys.stdout)    
        mqt.startInfo(self.parent.log)
#        textVBox = QVBoxLayout()
#        textVBox.addWidget(self.textBox)
        
        buttonHBox = QHBoxLayout()
        aboutButton = QPushButton("&ABOUT")
        self.connect(aboutButton, SIGNAL("clicked()"), self.on_about)
        buttonHBox.addWidget(aboutButton)
        saveButton = QPushButton("Save Printout")
        self.connect(saveButton, SIGNAL("clicked()"), self.on_save)
        buttonHBox.addWidget(saveButton)
        saveHTMLButton = QPushButton("Save HTML")
        self.connect(saveHTMLButton, SIGNAL("clicked()"), self.on_save_html)
        buttonHBox.addWidget(saveHTMLButton)
        
        layout = QVBoxLayout(self)
        layout.addWidget(self.textBox)
        layout.addLayout(buttonHBox)
        
    def on_save(self):
        printOutFilename, filt = QFileDialog.getSaveFileName(self,
                "Save as PRT file...", ".prt",
                "PRT files (*.prt)")

        self.textBox.selectAll()
        text = self.textBox.toPlainText()
        fout = open(printOutFilename,'w')
        fout.write(text)
        fout.close()

        self.textBox.append('Session saved to printout file:')
        self.textBox.append(printOutFilename)
        
    def on_save_html(self):
        self.parent.report.outputhtml()
        
    def on_about(self):
        dialog = mqt.AboutDlg(self)
        if dialog.exec_():
            pass
        
        
class PlotBlock(QWidget):
    def __init__(self, parent=None):
        super(PlotBlock, self).__init__(parent)
        self.parent = parent
        #self.plot_legend = False
        
        # Prepare plot window
        layout = QVBoxLayout(self)
        self.legendChB = QCheckBox("Show L&egend")
        self.legendChB.setChecked(True)
        #self.connect(self.legendChB, SIGNAL("stateChanged()"), self.on_something_changed)
        self.parent.canvas = mqt.MatPlotWin()
        canvastools = mqt.MatPlotTools(self.parent.canvas, self.parent)
        layout.addWidget(self.legendChB)
        layout.addWidget(self.parent.canvas)
        layout.addWidget(canvastools)
        
    def on_something_changed(self):
        self.plot_legend = False
        self.plot_legend = self.legendChB.isChecked()
        self.on_show()
        
#    def on_show(self, plotGuesses=False, 
#        plotFit=False, plotNorm=False, plotPooled=False):
#        self.parent.canvas.axes.clear()
#        self.parent.canvas.axes.grid(True)
#        
#            
#        if plotNorm:
#            for session in self.parent.fits:
#                self.parent.canvas.axes.semilogx(session.data.X, 
#                    session.data.normY, 'o', label=session.data.title)
#                logplotX = np.log10(session.data.X)
#                plotX = 10 ** np.linspace(np.floor(np.amin(logplotX) - 1),
#                    np.ceil(np.amax(logplotX)), 100)
#                plotYg = session.eq.equation(plotX, session.eq.normpars)
#                self.parent.canvas.axes.semilogx(plotX, plotYg, 'b-')
#                
#        elif plotPooled:
##            self.parent.canvas.axes.plot(self.parent.pooledfit.data.avX,
##                self.parent.pooledfit.data.avY, 'ro', label='average')
#            self.parent.canvas.axes.errorbar(self.parent.pooledfit.data.avX, 
#                self.parent.pooledfit.data.avY, 
#                yerr=self.parent.pooledfit.data.avS, fmt='o', ecolor='b', label='average')
#            logplotX = np.log10(self.parent.pooledfit.data.avX)
#            plotX = 10 ** np.linspace(np.floor(np.amin(logplotX) - 1),
#                np.ceil(np.amax(logplotX)), 100)
#            plotYg = self.parent.pooledfit.eq.equation(plotX,
#                self.parent.pooledfit.eq.pars)
#            self.parent.canvas.axes.semilogx(plotX, plotYg, 'b-')
#            
#        else:
#            for set in self.parent.data:
#                if set.S.any() == 0:
#                    self.parent.canvas.axes.semilogx(set.X, set.Y, 'o', label=set.title)
#                else: 
#                    self.parent.canvas.axes.errorbar(set.X, set.Y, yerr=set.S,
#                        fmt='o', label=set.title)
#                    self.parent.canvas.axes.set_xscale('log')
#            
#        if plotGuesses:
#            for session in self.parent.fits:
#                logplotX = np.log10(session.data.X)
#                plotX = 10 ** np.linspace(np.floor(np.amin(logplotX) - 1),
#                    np.ceil(np.amax(logplotX) + 1), 100)
#                plotYg = session.eq.equation(plotX, session.eq.pars)
#                self.parent.canvas.axes.semilogx(plotX, plotYg, 'y-')
#                
#        if plotFit:
#            for session in self.parent.fits:
#                logplotX = np.log10(session.data.X)
#                plotX = 10 ** np.linspace(np.floor(np.amin(logplotX) - 1),
#                    np.ceil(np.amax(logplotX) + 1), 100)
#                plotYg = session.eq.equation(plotX, session.eq.pars)
#                self.parent.canvas.axes.semilogx(plotX, plotYg, 'b-')
#                
#
#
#        if self.legendChB.isChecked():
#            self.parent.canvas.axes.legend(loc=2)
#            
#        self.parent.canvas.draw()
        
        
        
        
#
#        for row in range(self.series_list_model.rowCount()):
#            model_index = self.series_list_model.index(row, 0)
#            checked = self.series_list_model.data(model_index,
#                Qt.CheckStateRole) == Qt.Checked
#            name = str(self.series_list_model.data(model_index))
#
#            if checked:
#                has_series = True
#                self.data.add_series_to_fit(name)
#                set = self.data.get_series_data(name)
#                X = set[0]
#                Y = set[1]
#                Yerr = set[2]
#                self.axes.semilogx(X, Y, 'o', label=name)
#
#        if self.fitted:
#
#            self.equation.calcCurves()
#            xy = self.equation.curves
#            for i in range(self.equation.nsets):
#                #name1 = 'fit'+(i+1)
#                self.axes.semilogx(xy[0], xy[i+1], 'r-')
#
#        self.canvas.draw()
