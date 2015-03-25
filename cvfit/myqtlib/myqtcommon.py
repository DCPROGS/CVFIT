import sys
import socket
import datetime

try:
    from PyQt4.QtGui import *
    from PyQt4.QtCore import *
except:
    raise ImportError("pyqt module is missing")

try:
    from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
    from matplotlib.figure import Figure
    from matplotlib import scale as mscale
except:
    raise ImportError("matplotlib module is missing")

def startInfo(log):
    """
    Get date, time, machine info, etc.
    """
    log.write("CVFIT: Curve fitting program")
    log.write("Date and time of analysis: " + str(datetime.datetime.now())[:19])
    machine = socket.gethostname()
    system = sys.platform + '\n' + sys.version
    log.write("Machine: {0};   \nSystem: {1}".format(machine, system))
    

class PrintLog:
    """
    Write stdout to a QTextEdit.
    out1 = QTextEdit, QTextBrowser, etc.
    out2 = sys.stdout, file, etc.
    """
    def __init__(self, out1, out2=None):
        self.out1 = out1
        self.out2 = out2
    def write(self, text):
        self.out1.append(text.rstrip('\n'))
        if self.out2:
            self.out2.write(text)
                

class MatPlotWin(FigureCanvas):
    """
    """
    def __init__(self, size=(6.0, 4.0), fsize=8):
        # Prepare matplotlib plot window
        self.fig = Figure(size, dpi=85)
        self.axes = self.fig.add_subplot(111)
        self.axes.autoscale_view(True,True,True)
        self.fontsize = fsize
        for loc, spine in self.axes.spines.items(): #iteritems():
            if loc in ['right','top']:
                spine.set_color('none') # don't draw spine
        self.axes.xaxis.set_ticks_position('bottom')
        self.axes.yaxis.set_ticks_position('left')
        for label in self.axes.xaxis.get_ticklabels():
            label.set_fontsize(self.fontsize)
        for label in self.axes.yaxis.get_ticklabels():
            label.set_fontsize(self.fontsize)
#        self.mplTools = NavigationToolbar(self.canvas, self.parent)
        FigureCanvas.__init__(self, self.fig)        
        
class MatPlotTools(NavigationToolbar):
    """
    """
    def __init__(self, canvas, parent):
        NavigationToolbar.__init__(self, canvas, parent)
        
        
class AboutDlg(QDialog):
    def __init__(self, parent=None):
        QDialog.__init__(self)

        okButton = QPushButton("&OK")
        button_layout = QHBoxLayout()
        button_layout.addStretch()
        button_layout.addWidget(okButton)
        button_layout.addStretch()

        movie_screen = self.movie_screen()
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("<p align=center><b>Welcome to CVFIT: "
        "Curve fitting!</b></p>"))
        layout.addWidget(movie_screen)
        layout.addLayout(button_layout)

        self.connect(okButton, SIGNAL("clicked()"),
        self, SLOT("accept()"))
        self.setStyleSheet("QWidget { background-color: %s }"% "white")
        self.setWindowTitle("About CVFIT: Curve fitting")

    def movie_screen(self):
    # set up the gif movie screen
        movie_screen = QLabel()
        # expand and center the label
        movie_screen.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        movie_screen.setAlignment(Qt.AlignCenter)
        movie = QMovie("dca2.gif", QByteArray(), self)
        movie.setCacheMode(QMovie.CacheAll)
        movie.setSpeed(100)
        movie_screen.setMovie(movie)
        movie.start()
        return movie_screen