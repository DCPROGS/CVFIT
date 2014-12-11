import codecs

import numpy as np
import markdown


class Report(object):

    '''
    A class that makes markdown file and output it as html file
    '''

    def __init__(self, filename):
        self.filename = filename + '.md'
        self.file = open(self.filename, 'w')

    def title(self, titletext, titlenumber):
        # Add a header
        self.file.write('\n' + '#' * titlenumber + ' ' + titletext + '\n')

    def paragraph(self, paragraphtext):
        # Add a paragraph
        # The paragraphtext is a list which contains each line of string
        self.file.write('\n' + paragraphtext + '\n')

    def image(self, imagefile):
        self.file.write('\n' + '![Alt text](' + imagefile + ')' + '\n')

    def tabletitle(self, titlelist):
        writetable = np.repeat(titlelist, 2).astype('string')
        writetable[::2] = ' | '
        self.file.write('\n' + writetable.tostring() + ' | ' + '\n')
        writetable[1::2] = '-' * 20
        self.file.write(writetable.tostring() + ' | ' + '\n')

    def table(self, nparray):
        # Create table form numpy array
        if len(nparray.shape) == 2:
            for line in nparray:
                writetable = np.repeat(line, 2).astype('string')
                writetable[::2] = ' | '
                self.file.write(writetable.tostring() + ' | ' + '\n')

        else:
            writetable = np.repeat(nparray, 2).astype('string')
            writetable[::2] = ' | '
            self.file.write(writetable.tostring() + ' | ' + '\n')

    def outputhtml(self):
        self.file.close()
        input_file = codecs.open(self.filename, mode='r', encoding='utf-8')
        text = input_file.read()
        html = markdown.markdown(
            text, extensions=['markdown.extensions.tables'])

        output_file = codecs.open(self.filename[:-3] + '.html', 'w',
                                  encoding='utf-8',
                                  errors='xmlcharrefreplace')
        output_file.write(html)
        output_file.close()

    def close(self):
        self.file.close()



