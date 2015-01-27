import os
import sys
import time
import socket
import codecs
import platform

import numpy as np
import markdown


class Report(object):

    '''
    A class that makes markdown file and output it as html file
    '''

    def __init__(self, filename):
        
        self.path, self.filename = set_cwd('CVFIT_fits', filename)
        self.mdfile = 'report_' + self.filename + '.md'
        self.f = open(self.mdfile, 'w')
         
        self.title("CVFIT: Curve fitting program\n", 1)
        str = ("Date and time of analysis: %4d/%02d/%02d %02d:%02d:%02d\n"
            %time.localtime()[0:6])
        machine = socket.gethostname()
        system = sys.platform
        str += "Machine: %s; System: %s\n" %(machine, system)
        self.f.write(str)

    def title(self, titletext, titlenumber):
        # Add a header
        self.f.write('\n' + '#' * titlenumber + ' ' + titletext + '\n')

    def paragraph(self, paragraphtext):
        # Add a paragraph
        # The paragraphtext is a list which contains each line of string
        self.f.write('\n' + paragraphtext + '\n')

    def image(self, imagefile):
        self.f.write('\n' + '![Alt text](' + imagefile + ')' + '\n')

    def tabletitle(self, titlelist):
        writetable = np.repeat(titlelist, 2).astype('string')
        writetable[::2] = ' | '
        self.f.write('\n' + writetable.tostring() + ' | ' + '\n')
        writetable[1::2] = '-' * 20
        self.f.write(writetable.tostring() + ' | ' + '\n')

    def table(self, tabletext):
        # Create table form numpy array
        for line in tabletext:
            writetable = np.repeat(line, 2).astype('string')
            writetable[::2] = ' | '
            self.f.write(writetable.tostring() + ' | | ' + '\n')
           
    def dataset(self, name, set):
        self.paragraph(name)
        rows = set.split('\n')
        self.tabletitle(rows[1].split('\t'))
        tabletxt = []
        for i in range(2, len(rows)):
            tabletxt.append(rows[i].split('\t'))
        self.table(np.array(tabletxt))
            

    def outputhtml(self):
        self.f.close()
        input_file = codecs.open(self.mdfile, mode='r', encoding='utf-8')
        text = input_file.read()
        html = markdown.markdown(
            text, extensions=['markdown.extensions.nl2br', 'markdown.extensions.tables'])

        output_file = codecs.open(self.mdfile.split('.')[0] + '.html', 'w',
                                  encoding='utf-8',
                                  errors='xmlcharrefreplace')
        output_file.write(html)
        output_file.close()

    def close(self):
        self.f.close()



def set_cwd(dirname, filename):
    os.chdir(os.path.expanduser("~"))
    if platform.system() == 'Windows':
        os.chdir(filename.split('/')[0]+'//')
    if not os.path.isdir(dirname):
        os.mkdir(dirname)
    os.chdir(dirname)

    if not filename:
        filename = ("%4d%02d%02d_%02d%02d%02d" %time.localtime()[0:6])
    else:
        filename = filename.split('/')[-1].split('.')[0]
        
    timesuf = ("_%4d%02d%02d_%02d%02d%02d" %time.localtime()[0:6])
    if not os.path.isdir(filename+timesuf):
        os.mkdir(filename+timesuf)
    os.chdir(filename+timesuf)

    return os.getcwd(), filename