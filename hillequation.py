import sys
import os
import codecs


import numpy as np
from scipy import optimize, stats
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import markdown
import prettyplotlib as ppl


from toolkit import Report, check_input


class Cell(object):

    '''
    A class that store the inportant information for each cell
    '''

    def __init__(self, concentration, response):
        '''
        Make the class by specifying the concentration and response.
        Original data is then sorted according to the concentration.
        '''

        # Keep track of the original data
        self.originalconcentration = concentration
        self.originalresponse = response
        # Sort the data according to concentration
        result = np.vstack((concentration, response))
        result.sort()
        self.concentration = result[0]
        self.response = result[1]
        self.size = len(result[0])


class HillEquation(Cell):

    '''
    A class for fitting the data with Hill equation.
    '''

    def __init__(self, cell):
        Cell.__init__(self, cell.concentration, cell.response)

    def set_condition_hillequation(self):
        '''
        Setting the conditions for fitting.
        '''

    # Select weighting method
        self.weightmode = None
        for i in np.unique(self.concentration):
            if len(np.where(self.concentration == i)[0]) == 1:
                self.weightmode = 1
                print 'Concentration', i, 'has only one value.'
                print 'Using equal weighted method.'
                break

        if self.weightmode is None:
            print 'Please select the weight calculation option.'
            print 'Mode 1: Equally weighted.'
            print 'Mode 2: Weight by the standard error.(Default)'
            self.weightmode = check_input('Mode number:', ['1', '2'], 2)

        # Select the number of Components
        print 'Please key in the number of the Components. Default = 1'
        self.Component = check_input('Number of the Components:',
                                     ['1', '2'], 1)

        # Detect the trend
        # Selecting whether the trend is positive or negative
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            self.concentration, self.response)
        if slope > 0:
            trend = 'positive'
            self.trend = 1
        elif slope < 0:
            trend = 'negative'
            self.trend = -1
        print 'Select the slope.', '(default:', trend, ')'
        print 'Please key in 1 for positive slope and -1 for negative slope'
        self.trend = check_input('-1 or 1:', ['-1', '1'], self.trend)

        if self.trend == 1:
            # Select to fix Ymin/Y(0) or not
            print 'Please select to fix the Y(0) or not.'
            print '0 for False and 1 for True. (Default = 1)'
            self.YminFixed = check_input('0 or 1:', ['0', '1'], 1)
        else:
            self.YminFixed = 0

        # Set the initial value for fixing Ymax
        self.YmaxFixed = False
        return self.weightmode, self.Component, self.trend, \
            self.YminFixed, self.YmaxFixed

    def directset_condition_hillequation(self, weightmode, Component, trend,
                                         YminFixed, YmaxFixed):
        self.weightmode = weightmode
        self.Component = Component
        self.trend = trend
        self.YminFixed = YminFixed
        self.YmaxFixed = YmaxFixed

    def hill_equation_guess(self):
        '''
        Calculate the guess for fitting with hill equation.
        '''

        if self.Component == 1:
            self.guess = np.empty(4)
            # Four parameters will be needed to fit one Component Hill
            # equation

            if self.trend == 1:
                # If response increases with concentration

                # Determine Y(0)
                if self.YminFixed:
                    self.guess[0] = 0
                else:
                    self.guess[0] = np.mean(
                        self.response[self.concentration == self.concentration[0]])
                if self.YmaxFixed:
                    self.guess[1] = 1
                else:
                    # Determine Ymax
                    self.guess[1] = np.mean(
                        self.response[self.concentration == self.concentration[-1]]) \
                        - self.guess[0]
            else:
                # If response decreases with concentration

                # Determine Y(0)
                self.guess[0] = np.mean(
                    self.response[self.concentration == self.concentration[-1]])

                # Determine Ymin
                self.guess[1] = np.mean(self.response[self.concentration == self.concentration[0]]) - \
                    self.guess[0]

            # Determine Kr
            Kr = self.guess[2] = 10 ** ((np.log10(self.concentration[0]) +
                                         np.log10(self.concentration[-1])) / 2)

            # Determine nH
            Ymax = np.amax(self.response)
            LinRegressCon = self.concentration[self.response < Ymax]
            LinRegressResponse = self.response[self.response < Ymax]
            LinRegressX = np.log10(LinRegressCon) - np.log10(self.guess[2])
            LinRegressY = np.log10(
                (LinRegressResponse / Ymax) / (1 - (LinRegressResponse / Ymax)))
            slope, intercept, r_value, p_value, std_err = stats.linregress(
                LinRegressX, LinRegressY)
            nH = self.guess[3] = slope

        elif self.Component == 2:
            print 'Two Components fitting is not completed.'
            sys.exit(0)

    def fitting_hill_equation(self):
        '''
        Use least square fitting method to fit data with hill equation.
        '''

        # Calculate weight
        if self.weightmode == 1:
            self.weight = np.ones(self.size)
        elif self.weightmode == 2:
            self.weight = np.empty(self.size)
            for index, concentration in enumerate(cell.concentration):
                response = cell.response[
                    cell.concentration == concentration]
                weight[index] = 1 / \
                    (np.std(response) / np.sqrt(len(response)))

        # Least square fitting
        coeffs, params = optimize.leastsq(
            hill_equation_residuals, self.guess,
            args=(self.response, self.weight, self.concentration,
                  self.Component, self.YminFixed, self.YmaxFixed))
        self.coeffs = coeffs

    def normalise_hill_equation(self):
        '''
        Nomalise the coefficients and response.
        '''

        if self.trend == 1:
            # Nomalise the coefficients by fixing the Y(0) and Ymax
            self.normalised_coeffs = self.coeffs.copy()
            self.normalised_coeffs[0] = 0
            self.normalised_coeffs[1] = 1

            # Nomalise the response
            self.normalised_response = \
                (self.response.copy() -
                 self.coeffs[0]) / self.coeffs[1]
        elif self.trend == -1:
            # Nomalise the coefficients by fixing the Y(0) and Ymax
            self.normalised_coeffs = self.coeffs.copy()
            self.normalised_coeffs[0] = 1
            self.normalised_coeffs[1] = -1

            # Nomalise the response
            self.normalised_response = 1 - \
                (self.response.copy() -
                 self.coeffs[0]) / self.coeffs[1]

    def calculate_average_standarderror(self):
        '''
        Calculate the average response in one concentration
        and the stardard deviation of the mean
        '''

        self.uniqueconcentration = np.array([])
        self.averageresponse = np.array([])
        self.standarderror = np.array([])
        for con in np.unique(self.concentration):
            response = self.response[self.concentration == con]
            self.averageresponse = np.append(
                self.averageresponse, np.mean(response))
            self.uniqueconcentration = np.append(
                self.uniqueconcentration, con)
            self.standarderror = np.append(
                self.standarderror, np.std(response) / np.sqrt(len(response)))

    def plot_originaldata_guess_fit(self, filename, index):
        '''
        Plot the original data with guess and fitted curve
        in log graph.
        '''

        fig, ax = plt.subplots()
        # Plot original data
        ppl.plot(ax, self.concentration, self.response,
                 'o', label='experimental data')

        logplotX = np.log10(self.concentration)
        plotX = 10 ** np.linspace(np.floor(np.amin(logplotX) - 1),
                                  np.ceil(np.amax(logplotX) + 1), 100)
        self.xmin = plotX[0].copy()
        self.xmax = plotX[-1].copy()

        # Plot guess curve
        plotYguess = hill_equation(
            self.guess, plotX, self.Component, self.YminFixed, False)
        ppl.plot(ax, plotX, plotYguess, label='guess')

        # Plot fitted curve
        plotYfit = hill_equation(
            self.coeffs, plotX, self.Component, self.YminFixed, False)
        ppl.plot(ax, plotX, plotYfit, label='fit')

        ppl.legend(ax, loc='lower right')
        ax.set_xscale('log')
        ax.set_title(
            os.path.split(filename)[1][:-4] + ' ' + str(index + 1))
        fig.savefig(filename[:-4] + '_' + str(index + 1) + '.png')
        plt.close(fig)


def hill_equation(coeffs, concentration, Component, YminFixed, YmaxFixed):
    '''
    The hill equation
    '''

    if YminFixed:
        ymin = 0
    else:
        ymin = coeffs[0]

    if YmaxFixed:
        ymax = 1
    else:
        ymax = coeffs[1]

    equation = ymin + ymax * (
        ((concentration / coeffs[2]) ** coeffs[3])
        /
        (1 + (concentration / coeffs[2]) ** coeffs[3]))
    return equation


def hill_equation_residuals(coeffs, response, weight, concentration,
                            Component, YminFixed, YmaxFixed):
    '''
    Calculate the weighted residuals for the hill equation
    '''

    result = weight * \
        (response -
         hill_equation(coeffs, concentration,
                       Component, YminFixed, YmaxFixed))
    return result


def fitting_curve_hill_equation(filename, celllist, report):
    print 'Do you want to select fitting conditions separately?'
    print '0 means using one condition for all cells (Default)'
    print '1 means key in condition differently for each cell.'
    separate_condition = check_input('0 or 1:', ['0', '1'], 0)
    processlist = []
    for index, cell in enumerate(celllist):
        print 'Cell:', str(index + 1)
        # Creat HillEquation class
        cell = HillEquation(cell)
        processlist.append(cell)
        # Setting fitting conditions
        if (separate_condition == 0) and (index == 0):
            weightmode, Component, trend, YminFixed, YmaxFixed = cell.set_condition_hillequation()
        elif separate_condition == 0:
            cell.directset_condition_hillequation(
                weightmode, Component, trend, YminFixed, YmaxFixed)
        elif separate_condition == 1:
            cell.set_condition_hillequation()

        # Calculate the guess
        cell.hill_equation_guess()

        # Fitting the hill equation
        cell.fitting_hill_equation()

        # Normalise the result
        cell.normalise_hill_equation()

        # Plot the original data, guess and fit
        cell.plot_originaldata_guess_fit(filename, index)

        print 'Fitting completed'
        print 'Y(0)', cell.coeffs[0], 'Ymax', cell.coeffs[1],
        print 'Kr', cell.coeffs[2], 'nH', cell.coeffs[3]

    # Print all the coefficients into the report
    report.title('Fiting result:', 1)
    report.title('Cofficients:', 2)
    title = ['Y(0)', 'Ymax', 'Kr', 'nH']
    report.tabletitle(title)
    table = np.vstack([cell.coeffs for cell in processlist])
    report.table(table)

    # Plotting three graph
    # Generate the common x axis
    xmin = np.floor(np.log10(min([cell.xmin for cell in processlist])))
    xmax = np.ceil(np.log10(max([cell.xmax for cell in processlist])))
    plotx = 10 ** np.linspace(xmin, xmax, 100)
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    fig3, ax3 = plt.subplots()
    for index, cell in enumerate(processlist):
        # Plot original data and fitted curve
        ppl.plot(ax1, cell.concentration, cell.response,
                 'o', color=ppl.colors.set2[index])
        ppl.plot(ax1, plotx, hill_equation(cell.coeffs, plotx,
                                           cell.Component, YminFixed, False),
                 color=ppl.colors.set2[index], label=str(index + 1))

        # Plot normalised curve only
        ppl.plot(ax2, plotx, hill_equation(cell.normalised_coeffs, plotx,
                                           cell.Component, YminFixed, False),
                 color=ppl.colors.set2[index], label=str(index + 1))

        # Plot normalised curve with normalised data
        ppl.plot(ax3, cell.concentration, cell.normalised_response,
                 'o', color=ppl.colors.set2[index])
        ppl.plot(ax3, plotx, hill_equation(cell.normalised_coeffs, plotx,
                                           cell.Component, YminFixed, False),
                 color=ppl.colors.set2[index], label=str(index + 1))

    ppl.legend(ax1, loc='lower right')
    ax1.set_xscale('log')
    ax1.set_title(os.path.split(filename)[1][:-4])
    fig1.savefig(filename[:-4] + '_originaldata_fittedcurve' + '.png')
    plt.close(fig1)

    ppl.legend(ax2, loc='lower right')
    ax2.set_xscale('log')
    ax2.set_title(os.path.split(filename)[1][:-4])
    fig2.savefig(filename[:-4] + '_normalisedfittedcurve' + '.png')
    plt.close(fig2)

    ppl.legend(ax3, loc='lower right')
    ax3.set_xscale('log')
    ax3.set_title(os.path.split(filename)[1][:-4])
    fig3.savefig(
        filename[:-4] + '_normaliseddata_fittedcurve' + '.png')
    plt.close(fig3)

    # print three plots to report
    report.title('Fitted curve:', 2)
    htmlfilename = os.path.split(filename)[-1]
    report.title('Original data with fitted curve:', 3)
    report.image(
        htmlfilename[:-4] + '_originaldata_fittedcurve' + '.png')
    report.title('Normalised fitted curve:', 3)
    report.image(htmlfilename[:-4] + '_normalisedfittedcurve' + '.png')
    report.title('Normalised data and fitted curve:', 3)
    report.image(
        htmlfilename[:-4] + '_normaliseddata_fittedcurve' + '.png')

    # Pool all the data together and fit
    totalcell = HillEquation(Cell(
        np.hstack([cell.concentration for cell in processlist]),
        np.hstack([cell.normalised_response for cell in processlist])))
    totalcell.directset_condition_hillequation(
        weightmode, Component, trend, YminFixed, True)
    totalcell.hill_equation_guess()
    totalcell.fitting_hill_equation()
    totalcell.calculate_average_standarderror()

    # Plot the fiting of nomalised curve
    fig, ax = plt.subplots()
    ax.errorbar(totalcell.uniqueconcentration,
                totalcell.averageresponse,
                yerr=totalcell.standarderror,
                fmt='o', color=ppl.colors.set2[0],
                ecolor=ppl.colors.set2[1], label='Errorbar')
    ppl.plot(ax, plotx, hill_equation(totalcell.coeffs, plotx,
                                      totalcell.Component, YminFixed, True),

             color=ppl.colors.set2[2], label='Fitting')
    ppl.legend(ax, loc='lower right')
    ax.set_xscale('log')
    ax.set_title(os.path.split(filename)[1][:-4])
    fig.savefig(filename[:-4] + '_allfitting' + '.png')
    plt.close(fig)

    # print the final result
    report.title('Fitting with all normalised data:', 1)
    report.title('Fitting plot:', 2)
    report.image(htmlfilename[:-4] + '_allfitting' + '.png')
    report.title('Fitting cofficients:', 2)
    report.tabletitle(title)
    report.table(totalcell.coeffs)
    report.title('Fitting accuracy:', 2)
    title = ['Concentration', 'Average', 'Standard error']
    report.tabletitle(title)
    table = np.vstack((totalcell.uniqueconcentration,
                       totalcell.averageresponse,
                       totalcell.standarderror))
    table = np.transpose(table)
    report.table(table)

    report.outputhtml()
