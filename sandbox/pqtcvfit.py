#! /usr/bin/python

import time, os, sys, socket
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from array import array
#from ReadRandat import read_Data


__author__="remis"
__date__ ="$22-Jan-2010 22:36:24$"

class Form(QDialog):
    def __init__(self, parent=None):
        super(Form, self).__init__(parent)

        self.intro = "\n    CVFIT- general curve fitting program."
        self.inidict = {}
        self.data = []

        tab_widget = QTabWidget()
        tab1 = QWidget()
        tab1.setStyleSheet("QWidget { background-color: %s }"% "white")
        tab2 = QWidget()
        tab_widget.addTab(tab1, "Wellcome!")
        tab_widget.addTab(tab2, "CVFIT")

        ####### TAB 1 ##########
        movie_screen = self.movie_screen()
        tab1_layout = QVBoxLayout(tab1)
        tab1_layout.addWidget(QLabel("<p align=center><b>Welcome to DC_PyPs: "
        "Curve fitting!</b></p>"))
        tab1_layout.addWidget(movie_screen)
        tab1_layout.addWidget(QLabel("<p align=center><b>To start fitting go to "
        "CVFIT tab.</b></p>"))

        ####### TAB 2 ##########
        tab2_layout = QVBoxLayout(tab2)
        tab2_layout = self.cvfit_layout(tab2_layout)

        ##### Finalise main window ######
        vbox = QVBoxLayout()
        vbox.addWidget(tab_widget)
        quitButton = QPushButton("&QUIT")
        self.connect(quitButton, SIGNAL("clicked()"), self, SLOT("close()"))
        vbox.addWidget(quitButton)
        self.setLayout(vbox)
        self.setWindowTitle("DC_PyPs: Curve fitting")


    def start_info(self):
        """Display date, time, machine info, etc."""

        #self.str1 = "\npCVFIT.py is a clone of DC's CVFIT.FOR- general curve fitting program."
        #print self.str1
        str2 = "Date and time of analysis: %4d/%02d/%02d %02d:%02d:%02d" %time.localtime()[0:6]
        print str2 

        machine = socket.gethostname()
        system = sys.platform
        system_version = sys.version
        str3 = "Machine: %s; System: %s" %(machine, system)

        print "Machine: " + machine
        print "System: " + system + system_version

        return str2, str3


    #######   TAB 1: WELCOME!  START   ############
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
#######   TAB 1: WELCOME!  END   ############

#######   TAB 2: RANTEST FOR CONTINUOSLY VARIABLY DATA. START  #############
    def cvfit_layout(self, tab_layout):
        "Create Tab2 layout."
        tab_layout.addWidget(QLabel(self.intro))

        layout1 = QHBoxLayout()
        layout1.addWidget(QLabel("Configuration file:"))
        self.tb2e1 = QLineEdit()
        layout1.addWidget(self.tb2e1)
        self.tb2b1 = QPushButton("Browse for ini file...")
        layout1.addWidget(self.tb2b1)
        tab_layout.addLayout(layout1)

        layout2 = QHBoxLayout()
        layout2.addWidget(QLabel("Text data file:"))
        self.tb2e2 = QLineEdit()
        layout2.addWidget(self.tb2e2)
        self.tb2b2 = QPushButton("Browse for txt file...")
        layout2.addWidget(self.tb2b2)
        tab_layout.addLayout(layout2)

        self.tb2txt1 = QTextBrowser()
        str2, str3 = self.start_info()
        self.tb2txt1.append(str2)
        self.tb2txt1.append(str3)
        tab_layout.addWidget(self.tb2txt1)

        self.tb2b3 = QPushButton("Plot data")
        tab_layout.addWidget(self.tb2b3)

        self.tb2gr1 = QGraphicsView()
        scene = QGraphicsScene()
        scene.addLine(-400, 0, 270, 100)
        scene.addText("Hello")
        self.tb2gr1.setScene(scene)
        tab_layout.addWidget(self.tb2gr1)

        self.connect(self.tb2b1, SIGNAL("clicked()"), self.callback1)
        self.connect(self.tb2b2, SIGNAL("clicked()"), self.callback2)
        self.connect(self.tb2b3, SIGNAL("clicked()"), self.callback3)

        return tab_layout

    def callback3(self):
        pass

    def callback2(self):
        'Called by BROWSE FOR DATA TXT button in Tab2'
        #self.rancon_data, dfile = read_Data(self)
        file = QFileDialog.getOpenFileName(None,
         "Open Text Data File...", "", "TXT Files (*.txt)")
        self.tb2e2.setText(file)        
        data = self.read_txt_data(file)
        self.data = data
        self.tb2txt1.append('\nData loaded from: %s'%file)
        njset = len(data) / 3.0
        self.tb2txt1.append("File contains %d sets (X, Y, Yerr) of data"%njset)

    def read_txt_data(self, file):
        'Read data from a text file in which X and Y are in columns \
         and are Tab delimited.'
        f = open(file, 'r')
        lines =f.readlines()
        f.close()
        data1 = []
        for line in lines:
            line = line.strip("\n")    #remove newline
            values=line.split('\t')    #divide lines into values at tabs
            data1.append(values)
        print data1
        data = []
        for i in range(0, len(data1[0])):
            col = []
            for j in range(0, len(data1)):
                try:
                    col.append(float(data1[j][i]))
                except:
                    print 'no value'
            data.append(col)
        print data
        return data

    def callback1(self):
        'Called by BROWSE FOR INI button in Tab2'

        file = QFileDialog.getOpenFileName(None,
         "Open Configuration File...", "", "INI Files (*.ini *.INI)")
        self.tb2e1.setText(file)
        try:
            self.read_ini(file, 1)
            self.tb2txt1.append('\nIni file loaded: %(inifile)s \
            \n ndev = %(ndev)s \n ifile1 = %(ifile1)d \n ifitmode = %(ifitmode)d \
            \n ilog = %(ilog)d \n idiskq = %(idiskq)d \n titlef = %(titlef)s \
            \n infil = %(infil)s \n niobs %(niobs)d \n njset %(njset)d \
            \n ascfil = %(ascfil)s \n nsfit = %(nsfit)d'%self.inidict)
        except:
            self.tb2txt1.append("\nFailed to load INI file.")

    def read_ini(self, file, verbose=0):
        """Chose and read ini file."""
        self.inidict = {}
        self.inidict['inifile'] = file
        if verbose: print "will read mech file:", file

        ints = array('i')
        f = open(file, 'r')

        self.inidict['ndev'] = f.read(2)
        if verbose: print "ndev=", self.inidict['ndev']
        ints.fromfile(f,1)
        self.inidict['ifile1'] = ints.pop()
        if verbose: print "ifile1=", self.inidict['ifile1']
        ints.fromfile(f,1)
        self.inidict['ifitmode'] = ints.pop()
        if verbose: print "ifitmode=", self.inidict['ifitmode']
        ints.fromfile(f,1)
        self.inidict['ilog'] = ints.pop()
        if verbose: print "ilog=", self.inidict['ilog']
        ints.fromfile(f,1)
        self.inidict['idiskq'] = ints.pop()
        if verbose: print "idiskq=", self.inidict['idiskq']
	self.inidict['titlef'] = f.read(60)
        if verbose: print "titlef=", self.inidict['titlef']
        self.inidict['infil'] = f.read(33)
        if verbose: print "infil=", self.inidict['infil']
        #niobs- Maximimum number of observations per set
        ints.fromfile(f,1)
        self.inidict['niobs'] = ints.pop()
        if verbose: print "niobs=", self.inidict['niobs']
        #njset- Maximimum array size for sets (# of sets +5)?
        ints.fromfile(f,1)
        self.inidict['njset'] = ints.pop()
        if verbose: print "njset=", self.inidict['njset']
	self.inidict['ascfil'] = f.read(33)
        if verbose: print "ascfil=", self.inidict['ascfil']
        ints.fromfile(f,1)
        self.inidict['nsfit'] = ints.pop()
        if verbose: print "nsfit=", self.inidict['nsfit']

        f.close()
        return self.inidict

        

#######   TAB 2: RANTEST FOR CONTINUOSLY VARIABLY DATA. START  #############

    

if __name__ == "__main__":

    app = QApplication(sys.argv)
    form = Form()
    form.show()
    app.exec_()
