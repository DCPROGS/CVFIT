import sys
import os
import codecs
import multiprocessing

import numpy as np
from scipy import optimize, stats
import matplotlib.pyplot as plt
import markdown

import prettyplotlib as ppl


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


class HillEquation(object):

    '''
    A class that store the inportant information for each cell
    '''

    def __init__(self, concentration, response):
        '''
        Make the class by specifying the concentration and response
        '''

        # Keep track of the original data
        self.originalconcentration = concentration
        self.originalresponse = response
        # Sort the data according to concentration
        result = np.vstack((concentration, response))
        result.sort()
        self.concentration = result[0]
        self.response = result[1]

    def set_coeffs(self, coeffs):
        '''
        Setting the coefficient from the fitting process
        '''

        self.coeffs = coeffs
        # Nomalise the coefficients by fixing the Y(0) and Ymax
        normalised = self.coeffs.copy()
        normalised[0] = 0
        normalised[1] = 1
        self.normalised_coeffs = normalised
        # Nomalise the response
        self.normalised_response = \
            (self.response.copy() -
             self.coeffs[0]) / self.coeffs[1]

    def set_conres(self, concentration, response):
        '''
        Resetting the concentration and response
        '''

        self.concentration = concentration
        self.response = response

    def cal_avgstd(self):
        '''
        Calculate the average response in one concentration
        and the stardard deviation of the mean
        '''

        self.average = np.array([])
        self.avgcon = np.array([])
        self.constd = np.array([])
        for con in np.unique(self.concentration):
            response = self.response[self.concentration == con]
            self.average = np.append(self.average, np.mean(response))
            self.avgcon = np.append(self.avgcon, con)
            self.constd = np.append(
                self.constd, np.std(response) / np.sqrt(len(response)))


def check_input(text, accept, default):
    '''
    Check if the input is in acceptable range or not
    If not, ask to key in another value
    '''

    inputnumber = raw_input(text)
    if inputnumber:
        while inputnumber not in accept:
            print text
            inputnumber = raw_input(text)
    else:
        inputnumber = default

    return int(inputnumber)


def hill_equation_one(con, coeffs, YminFixed, YmaxFixed):
    '''
    Hill equation with only one component
    coeffs[0] = Y(0)
    coeffs[1] = Ymax
    coeffs[2] = Kr
    coeffs[3] = nH
    When YminFixed is True, Ymin is fixed to one
    When YmaxFixed is True, Ymax is fixed to one
    '''

    equation = coeffs[0] * (not YminFixed) + \
        (coeffs[1] * (not YmaxFixed) + YmaxFixed) * \
        ((con / coeffs[2]) ** coeffs[3]) \
        / \
        (1 + (con / coeffs[2]) ** coeffs[3])
    return equation


def hill_equation_two(con, coeffs, YminFixed, YmaxFixed):
    '''
    Hill equation with two component(incomplete)
    coeffs[0] = Y(0)
    coeffs[1] = Ymax
    coeffs[2] = Kr
    coeffs[3] = nH
    When YminFixed is True, Ymin is fixed to one
    When YmaxFixed is True, Ymax is fixed to one
    '''

    equation = coeffs[0] * (not YminFixed) + \
        (coeffs[1] * (not YmaxFixed) + YmaxFixed) * (
        ((con / coeffs[2]) ** coeffs[3])
        /
        (1 + (con / coeffs[2]) ** coeffs[3])
        +
        ((con / coeffs[4]) ** coeffs[5])
        /
        (1 + (con / coeffs[4]) ** coeffs[5]))
    return equation


def residuals(coeffs, response, con, weight, component, YminFixed, YmaxFixed):
    '''
    Calculate the weighted residuals for the hill equation
    '''

    if component == 1:
        result = weight * \
            (response -
             hill_equation_one(con, coeffs, YminFixed, YmaxFixed))
    elif component == 2:
        print 'Two components fitting is not completed'
        sys.exit(1)
        result = weight * \
            (response -
             hill_equation_two(con, coeffs, YminFixed, YmaxFixed))
    return result


def leastsqfit_hill(cell, component, YminFixed, weightmode, YmaxFixed):
    '''
    Use weighted least square method to fit the data to the hill equation
    '''

    # Calculate the weight
    if weightmode == 1:
        weight = np.ones(len(cell.concentration))
    elif weightmode == 2:
        weight = np.empty(len(cell.concentration))
        for index, concentration in enumerate(cell.concentration):
            response = cell.response[cell.concentration == concentration]
            weight[index] = 1 / \
                (np.std(response) / np.sqrt(len(response)))

    # Calculte the guess
    Ymin = 0
    if YmaxFixed:
        Ymax = 1
    else:
        Ymax = np.max(cell.response)

    Kr = 10 ** ((np.log10(cell.concentration[0]) +
                 np.log10(cell.concentration[-1])) / 2)

    LinRegressCon = cell.concentration[cell.response < Ymax]
    LinRegressResponse = cell.response[cell.response < Ymax]
    LinRegressX = np.log10(LinRegressCon) - np.log10(Kr)
    LinRegressY = np.log10(
        (LinRegressResponse / Ymax) / (1 - (LinRegressResponse / Ymax)))
    slope, intercept, r_value, p_value, std_err = stats.linregress(
        LinRegressX, LinRegressY)
    nH = slope

    # Fit the data
    guess = [Ymin, Ymax, Kr, nH]
    coeffs, params = optimize.leastsq(
        residuals, guess,
        args=(cell.response, cell.concentration,
              weight, component, YminFixed, YmaxFixed))
    cell.set_coeffs(coeffs)
    return guess, coeffs


def showplot(cell, index, component, YminFixed, weightmode):
    '''
    The aim of the function is to provide space for other evaluating process
    The current function is drawing a graph showing the fitting and the guess
    '''

    guess, coeffs = leastsqfit_hill(
        cell, component, YminFixed, weightmode, False)

    fig, ax = plt.subplots()
    ppl.plot(ax, cell.concentration, cell.response,
             'o', label='experimental data')
    logplotX = np.log10(cell.concentration)
    plotX = 10 ** np.linspace(np.floor(np.amin(logplotX) - 1),
                              np.ceil(np.amax(logplotX) + 1), 100)
    xmin = plotX[0].copy()
    xmax = plotX[-1].copy()
    plotYguess = hill_equation_one(plotX, guess, YminFixed, False)
    ppl.plot(ax, plotX, plotYguess, label='guess')
    plotYfit = hill_equation_one(plotX, coeffs, YminFixed, False)
    ppl.plot(ax, plotX, plotYfit, label='fit')
    ppl.legend(ax, loc='lower right')
    ax.set_xscale('log')
    ax.set_title(os.path.split(filename)[1][:-4] + ' ' + str(index + 1))
    fig.savefig(filename[:-4] + '_' + str(index + 1) + '.png')
    plt.close()
    print 'cell', index + 1, 'fitting complete.', \
        'Kr =', coeffs[2], 'nH =', coeffs[3]
    return coeffs, xmin, xmax


def valleyplot(totalcell, processlist):
    '''
    under testing
    '''

    nH, Kr = (np.array([]) for i in range(2))
    for index, cell in enumerate(processlist):
        Kr = np.append(Kr, cell.coeffs[2])
        nH = np.append(nH, cell.coeffs[3])
    nHmin = np.amin(nH)
    nHmax = np.amax(nH)
    Krmin = np.amin(Kr)
    Krmax = np.amax(Kr)

    Krrange = np.linspace(Krmin, Krmax, 100)
    nHrange = np.expand_dims(np.linspace(nHmin, nHmax, 100), axis=1)
    response = np.zeros((100, 100))
    for response, concentration in zip(totalcell.response, totalcell.concentration):
        response += (response - \
            ((concentration / Krrange) ** nHrange) \
            / \
            (1 + (concentration / Krrange) ** nHrange))**2
    fig, ax = plt.subplots()
    ppl.pcolormesh(fig, ax, response)
    ax.set_title(os.path.split(filename)[1][:-4])
    fig.savefig(filename[:-4] + '_pcolor' + '.png')
    plt.close()


print 'Please type in the loaction of the csv file'
filename = raw_input('filename:').strip()
while (os.path.exists(filename) is False) or (filename[-4:] != '.csv'):
    if os.path.exists(filename) is False:
        print filename, 'not found. Please type in again.'
    elif filename[-4:] != '.csv':
        print filename, 'is not a CSV file. Please type in again.'
    filename = raw_input('filename:').strip()

rawresult = np.genfromtxt(filename, delimiter=',')
setnum = rawresult.shape[1] / 2
print str(setnum), 'sets found.'

print 'Please type in the mode number:'
print 'Mode 1: Do different fitting for different data sets.(Default)'
print 'Mode 2: Do one fitting for all the data sets.'
rearange = check_input('Mode number:', ['1', '2'], 1)

if rearange == 2:
    result = rawresult[np.isfinite(rawresult)]
    processlist = [HillEquation(result[::2], result[1::2])]
elif rearange == 1:
    processlist = []
    for index in range(setnum):
        real = np.isfinite(rawresult[:, index * 2])
        cell = HillEquation(rawresult[real, index * 2],
                            rawresult[real, index * 2 + 1])
        processlist.append(cell)

weightmode = 2
breaker = None
for check in processlist:
    for i in np.unique(check.concentration):
        if len(np.where(check.concentration == i)[0]) == 1:
            weightmode = 1
            print 'Concentration', i, 'has only one value.'
            print 'Using equal weighted method.'
            breaker = True
            break
    if breaker:
        break

if weightmode == 2:
    print 'Please select the weight calculation option.'
    print 'Mode 1: Equally weighted.'
    print 'Mode 2: Weight by the standard error.(Default)'
    weightmode = check_input('Mode number:', ['1', '2'], 2)

print 'Please key in the number of the components. Default = 1'
component = check_input('Number of the components:', ['1', '2'], 1)

print 'Please select to fix the Y(0) or not.'
print '0 for False and 1 for True. (Default = 1)'
YminFixed = check_input('0 or 1:', ['0', '1'], 1)


sumfit = np.empty([setnum, 4])
sumxmin = np.empty(setnum)
sumxmax = np.empty(setnum)
for index, cell in enumerate(processlist):
    coeffs, xmin, xmax = showplot(
        cell, index, component, YminFixed, weightmode)
    sumxmin[index] = xmin
    sumxmax[index] = xmax

xmin = np.floor(np.log10(np.amin(sumxmin)))
xmax = np.ceil(np.log10(np.amax(sumxmax)))
plotx = 10 ** np.linspace(xmin, xmax, 100)


totalcell = HillEquation(np.array([]), np.array([]))
concentration = np.array([])
normalised_response = np.array([])
for cell in processlist:
    concentration = np.append(concentration, cell.concentration)
    normalised_response = np.append(
        normalised_response, cell.normalised_response)

totalcell.set_conres(concentration, normalised_response)
guess, coeffs = leastsqfit_hill(
    totalcell, component, YminFixed, weightmode, True)
totalcell.cal_avgstd()
fig, ax = plt.subplots()
ax.errorbar(totalcell.avgcon, totalcell.average,
            yerr=totalcell.constd, fmt='o', color=ppl.colors.set2[0],
            ecolor=ppl.colors.set2[1], label='Errorbar')
ppl.plot(ax, plotx, hill_equation_one(plotx, totalcell.coeffs,
                                      YminFixed, False),
         color=ppl.colors.set2[2], label='Fitting')
ppl.legend(ax, loc='lower right')
ax.set_xscale('log')
ax.set_title(os.path.split(filename)[1][:-4])
fig.savefig(filename[:-4] + '_fourth' + '.png')


fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()
fig3, ax3 = plt.subplots()
for index, cell in enumerate(processlist):
    ppl.plot(ax1, cell.concentration, cell.response,
             'o', color=ppl.colors.set2[index])
    ppl.plot(ax1, plotx, hill_equation_one(plotx, cell.coeffs,
                                           YminFixed, False),
             color=ppl.colors.set2[index], label=str(index + 1))

    ppl.plot(ax2, plotx, hill_equation_one(plotx, cell.normalised_coeffs,
                                           YminFixed, False),
             color=ppl.colors.set2[index], label=str(index + 1))

    ppl.plot(ax3, cell.concentration, cell.normalised_response,
             'o', color=ppl.colors.set2[index])
    ppl.plot(ax3, plotx, hill_equation_one(plotx, cell.normalised_coeffs,
                                           YminFixed, False),
             color=ppl.colors.set2[index], label=str(index + 1))

ppl.legend(ax1, loc='lower right')
ax1.set_xscale('log')
ax1.set_title(os.path.split(filename)[1][:-4])
fig1.savefig(filename[:-4] + '_first' + '.png')

ppl.legend(ax2, loc='lower right')
ax2.set_xscale('log')
ax2.set_title(os.path.split(filename)[1][:-4])
fig2.savefig(filename[:-4] + '_second' + '.png')

ppl.legend(ax3, loc='lower right')
ax3.set_xscale('log')
ax3.set_title(os.path.split(filename)[1][:-4])
fig3.savefig(filename[:-4] + '_third' + '.png')

plt.close()

valleyplot(totalcell, processlist)


report = Report(filename[:-4])
report.title('Original data:', 1)
report.paragraph('Number of cells been examined: ' + str(setnum))
report.title('Concentration and response:', 2)
title = np.empty(setnum * 2, dtype='S20')
title[::2] = 'Concentration:'
title[1::2] = 'Response:'
report.tabletitle(title)
rawresult = np.genfromtxt(filename, delimiter=',')
mask = np.isnan(rawresult)
rawresult = rawresult.astype('string')
rawresult[mask] = ' '
report.table(rawresult)
report.title('Fitted curve:', 1)
htmlfilename = os.path.split(filename)[-1]
report.title('Original data with fitted curve:', 2)
report.image(htmlfilename[:-4] + '_first' + '.png')
report.title('Normalised fitted curve:', 2)
report.image(htmlfilename[:-4] + '_second' + '.png')
report.title('Normalised data and fitted curve:', 2)
report.image(htmlfilename[:-4] + '_third' + '.png')
report.title('Cofficients:', 2)
title = ['Y(0)', 'Ymax', 'Kr', 'nH']
report.tabletitle(title)
table = np.empty([setnum, 4])
for index, cell in enumerate(processlist):
    table[index, :] = cell.coeffs
report.table(table)
report.title('Fitting with all normalised data:', 1)
report.title('Fitting plot:', 2)
report.image(htmlfilename[:-4] + '_fourth' + '.png')
report.title('Fitting cofficients:', 2)
report.tabletitle(title)
report.table(totalcell.coeffs)
report.title('Fitting accuracy:', 2)
title = ['Concentration', 'Average', 'Standard error']
report.tabletitle(title)
table = np.vstack((totalcell.avgcon, totalcell.average, totalcell.constd))
table = np.transpose(table)
report.table(table)


report.outputhtml()
