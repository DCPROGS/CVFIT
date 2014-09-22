CVFIT
=====

A curve fitting program which fits the data with Hill
equation with weighted least square and obtaining the coefficients
including Y(0), Ymax, nH and Kr.

#Features
This program can fits different cells with different hill equation
and generates the fitting result.
![Alt text](Example/Example_first.png)
All the coefficients can be normalised and produces several curves
with the same Ymax and Y(0).
![Alt text](Example/Example_second.png)
All the original data can also be normalised and generate a curve
which takes into account all the data.
![Alt text](Example/Example_third.png)
The standard deviation of the error in each concentration can be
plotted as error bar.
![Alt text](Example/Example_fourth.png)

#Requirements
**numpy, scipy, matplotlib** Widely used library when dealling with
scientific calculation

**Markdown** A library which convert Markdown file to html

    pip install markdown

**Prettyplotlib** A matplotlib-enhancer library

    pip install prettyplotlib
