import os
import xlrd
import numpy as np

from PyQt4.QtGui import *
from PyQt4.QtCore import *
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt

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
        
    def update_data(self):
        self.parent.data = []
        for row in range(self.dataListModel.rowCount()):
            i = self.dataListModel.index(row, 0)
            checked = (self.dataListModel.data(i, Qt.CheckStateRole) ==
                Qt.Checked)
            if checked:
                self.parent.data.append(self.allsets[row])
        
    def on_plot(self):
        self.update_data()
        #self.parent.plotblk.on_show()
#        plots.plot(self.parent.data, axes=self.parent.canvas.axes,
#            legend=self.parent.plotblk.legendChB.isChecked())
        plots.plot(self.parent.data, fig=self.parent.figure, 
            logX=self.parent.plotblk.logXChB.isChecked(),
            logY=self.parent.plotblk.logYChB.isChecked(),
            legend=self.parent.plotblk.legendChB.isChecked())
        self.parent.canvas.draw()

    def on_load(self):
        filename = QFileDialog.getOpenFileName(self.parent,
                "Open a data file...", ".",
                "Excel files (*.xls *.xlsx);;CSV files (*.csv);;TXT files (*.txt);;All files (*.*)")
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
        '(1) Hill equation',
        '(2) Langmuir equation',
        '(3) Straight line',
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
        
        self.parent.fits.pool(norm=True, output=self.parent.log)
        self.parent.fits.pooled.fit()
        self.parent.fits.pooled.calculate_errors()
        self.parent.fits.pooled.data.average_pooled()

        self.parent.log.write('\n***********************\n**************************')
        self.parent.log.write('\tNormalised and pooled data fit finished')
        self.parent.log.write(self.parent.fits.pooled.string_estimates())
        self.parent.log.write(self.parent.fits.pooled.string_liklimits())
        
        fplots = self.parent.fits.prepare_fplot('pooled')
        plots.plot([self.parent.fits.pooled.data], 
            fig=self.parent.figure, fplotsets=fplots, fplotline='b-',
            logX=self.parent.plotblk.logXChB.isChecked(),
            logY=self.parent.plotblk.logYChB.isChecked(),
            legend=self.parent.plotblk.legendChB.isChecked(), pooled=True)
        self.parent.canvas.draw()
        fname = 'fitted_normalised_pooled.png'
        self.parent.figure.savefig(fname)

        self.parent.report.title('Pooled data fit finished', 1)
        self.parent.report.paragraph(self.parent.fits.pooled.string_estimates())
        self.parent.report.paragraph(self.parent.fits.pooled.string_liklimits())
        self.parent.report.image(fname)
        
        fout = open(fname[:-4] + '1.txt', 'w')
        for i in range(len(self.parent.fits.pooled.data.avX)):
            fout.write('{0:.6e}\t{1:.6e}\t{2:.6e}\n'.
                format(self.parent.fits.pooled.data.avX[i], 
                self.parent.fits.pooled.data.avY[i],
                self.parent.fits.pooled.data.avS[i]))
        fout.close()
        fout = open(fname[:-4] + '2.txt', 'w')
        for i in range(len(fplots[0][0])):
            fout.write('{0:.6e}\t{1:.6e}\n'.
                format(fplots[0][0][i], fplots[0][1][i]))
        fout.close()

        
    def on_normalise(self):
        for session in self.parent.fits.list:
            session.eq.normalise(session.data)

        fname = 'all_fitted_normalised_curves.png'
        fplots = self.parent.fits.prepare_fplot('norm')
        plots.plot(self.parent.data, fig=self.parent.figure, 
            fplotsets=fplots, fplotline='b-',
            logX=self.parent.plotblk.logXChB.isChecked(),
            logY=self.parent.plotblk.logYChB.isChecked(),
            legend=self.parent.plotblk.legendChB.isChecked(), norm=True)
        self.parent.canvas.draw()
        self.parent.figure.savefig(fname)
        self.parent.report.title('Data normalised to the fitted maxima', 1)
        self.parent.report.image(fname)
        
    def on_fit(self):

        self.parent.report = Report(self.parent.fname)
        self.parent.report.title('Original data:', 1)
        self.parent.report.paragraph('Number of datasets loaded: ' + str(len(self.parent.fits.list)))
        for fit in self.parent.fits.list:
            self.parent.report.dataset(fit.data.title, str(fit.data))

        progressDlg = QProgressDialog('Fitting data set {0:d}'.format(1),
                                 "Cancel", 1, len(self.parent.fits.list))
        progressDlg.setWindowTitle('Fitting...')
        i = 1
        for fs in self.parent.fits.list:
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
                
        self.parent.log.write('\n\nAverage of all fits:')
        self.parent.log.write(self.parent.fits.string_average_estimates())
        self.parent.report.title('Average of all fits:', 1)
        self.parent.report.paragraph(self.parent.fits.string_average_estimates())

        fname = 'all_fittedcurves.png'
        fplots = self.parent.fits.prepare_fplot('fit')
        plots.plot(self.parent.data, fig=self.parent.figure, 
            fplotsets=fplots, fplotline='b-',
            logX=self.parent.plotblk.logXChB.isChecked(),
            logY=self.parent.plotblk.logYChB.isChecked(),
            legend=self.parent.plotblk.legendChB.isChecked())
        self.parent.canvas.draw()
        
#        plots.plot(self.parent.data, self.parent.fits, axes=self.parent.canvas.axes, plotFit=True, 
#            legend=self.parent.plotblk.legendChB.isChecked()) #, save_fig_name=fname)
#        self.parent.canvas.draw()
        self.parent.figure.savefig(fname)
        self.parent.report.image(fname)
        
    def plot_guess(self):
        fplots = self.parent.fits.prepare_fplot('guess')
        plots.plot(self.parent.data, fig=self.parent.figure, 
            fplotsets=fplots, fplotline='y-',
            logX=self.parent.plotblk.logXChB.isChecked(),
            logY=self.parent.plotblk.logYChB.isChecked(),
            legend=self.parent.plotblk.legendChB.isChecked())
        self.parent.canvas.draw()
        
    def on_guess(self):
        if self.parent.fits:
            self.plot_guess()
            dialog = dialogs.GuessDlg(self.parent.fits)
            if dialog.exec_():
                self.parent.fits = dialog.return_guesses()
            self.plot_guess()
        else:
            dialog = dialogs.WarningDlg('Please, load equation first!')
            if dialog.exec_():
                pass
        
    def on_equation(self):
        self.parent.datablk.update_data()
        self.parent.log.write('**********************************************')
        row = self.eqList.currentRow()
        if row == 0 or row == -1:
            self.parent.eqtype = 'Hill'
        elif row == 1:
            self.parent.eqtype = 'Langmuir'
        elif row == 2:
            self.parent.eqtype = 'Linear'
        else:
            self.parent.log.write("This eqation is not implemented yet.\n" +
                "Please, choose other equation.")
        dialog = dialogs.EquationDlg(self.parent.data, self.parent.eqtype, 
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
        layout = QVBoxLayout(self)
        
        checkHBox = QHBoxLayout()
        self.legendChB = QCheckBox("Show L&egend")
        self.legendChB.setChecked(True)
        #self.connect(self.legendChB, SIGNAL("stateChanged()"), self.on_something_changed)
        checkHBox.addWidget(self.legendChB)
        self.logXChB = QCheckBox("LogX")
        self.logXChB.setChecked(False)
        checkHBox.addWidget(self.logXChB)
        self.logYChB = QCheckBox("LogY")
        self.logYChB.setChecked(False)
        checkHBox.addWidget(self.logYChB)
        layout.addLayout(checkHBox)
        
        self.parent.figure = plt.figure()
        self.parent.canvas = FigureCanvas(self.parent.figure)
        #self.parent.canvas = mqt.MatPlotWin()
        canvastools = mqt.MatPlotTools(self.parent.canvas, self.parent)
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
