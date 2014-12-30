#! /usr/bin/python

__author__="remis"
__date__ ="$23-Feb-2010 10:20:30$"

import numpy as np

class DataSeries(object):
    """ A wrapper over a dictionary that holds data sets.
        Each series has a name and a list of X, Y and Yerr as its data.
        The length of all series is assumed to be the same?
        The sets can be read from CSV, TXT or DAT file.
        CSV:
        TXT (tab delimited as saved in Excel): each series contains three
        columns- X, Y, Yerr. First column first line- set title.
        DAT:
    """

    def __init__(self, filename=None):
        #self.load_from_txt_file(filename)
        self.data = {}
        self.names = []
        self.max_len = 0
        self.series_to_fit = []
        self.changed = 0
        self.nsets = len(self.series_to_fit)




    def clear_series_to_fit(self):
        ""
        self.series_to_fit = []

    def add_series_to_fit(self, name):
        ""
        self.series_to_fit.append(name)
        self.nsets = len(self.series_to_fit)

    def load_from_txt_file(self, filename=None):
        'Read data from a text file in which X and Y are in columns \
         and are Tab delimited.'

        f = open(filename, 'r')
        firstline = f.readline().strip("\n")
        names = firstline.split("\t")
        for i in range(int(len(names)/3.0)):
            self.names.append(names[3*i])

        lines =f.readlines()
        f.close()
        data1 = []
        for line in lines:
            line = line.strip("\n")    #remove newline
            values=line.split('\t')    #divide lines into values at tabs
            data1.append(values)
        #print data1
        data = []
        for i in range(0, len(data1[0])):
            col = []
            for j in range(0, len(data1)):
                try:
                    col.append(float(data1[j][i]))
                except:
                    #print 'no value'
                    pass
            data.append(col)

        sets = int(len(data) / 3.0)
        for i in range(sets):
            data0 = np.zeros((len(data[3*i])))
            dataout = np.array([data[3*i], data[3*i+1],
                                data[3*i+2], data0])

            
            #, data[3*i+1], data[3*i+2], np.zeros[len(data[3*i])]])

            #print 'dataout', dataout
#            dict = {}
#            dict['x'] = data[3*i]
#            dict['y'] = data[3*i+1]
#            dict['yerr'] = data[3*i+2]
#            dict['weights'] = np.zeros[len(data[3*i])]
            self.data[names[3*i]] = dataout
            if len(data[3*i]) > self.max_len: self.max_len = len(data[3*i])

    def load_from_csv_file(self, filename=None):

        if filename:
            for line in csv.reader(open(filename, 'rb')):
                self.names.append(line[0])
                self.data[line[0]] = map(int, line[1:])
                #self.datalen = len(line[1:])

    def series_names(self):
        """ Names of the data series
        """
        return self.names

#    def series_len(self):
#        """ Length of a data series
#        """
#        return self.datalen

    def series_count(self):
        return len(self.data)

    def get_series_data(self, name):
        return self.data[name]


if __name__ == "__main__":
    print "Hello World";
