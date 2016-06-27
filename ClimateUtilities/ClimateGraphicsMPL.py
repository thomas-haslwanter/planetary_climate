'''
Graphics interface customized for MatPlotLib

This is imported into ClimateUtilities. 
plotObj : class
          Mostly for saving plots in different formats
          
plot : function
       Takes a Curve object (defined in "ClimateUtilities" as input
       and produces a line plot. returns a "plotObj" object to allow
       further manipulation.

contour : function
        produces contour plots of an array. 

Change Log
----------
1/25/2012: Fixed bug in legends that made problem with Enthought 7.2
11/12/2013: plot()--Implemented plotting of multiple curves with
                differing resolution
06/15/2016: cleaned up, and ported to Python 3 [ThH]                

License
-------
BSD 3-clause (see https://www.w3.org/Consortium/Legal/2008/03-bsd-license.html)
'''

import numpy as np

#Try to import matplotlib plotting routines 
try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    
    #Configure the backend so things work properly
    #with the Python intepreter and idle. Not needed if you are
    #using ipython. In IPython you can eliminate
    #the following two commands.
    mpl.rcParams['backend'] = 'TkAgg'
    mpl.rcParams['interactive'] = True
    
except:
    print('matplotlib not found on your system')
    print('You can still run the courseware, but')
    print('will not be able to plot from inside Python')

class plotObj:
    '''
    A little class for use as a return object from plot(), so
    the user has an easy way to delete a window or save a plot
    to an eps file. 
    '''
    
    def __init__(self, workstation, plot, WorkstationResources=None):
        self.workstation = workstation
        self.plot = plot
        self.WorkstationResources = WorkstationResources
        
    def delete(self):
        '''Deletes a plot window '''
        print("In MatPlotLib just click the goaway button to dispose a plot")
        
    def save(self, filename='plot.png'):
        '''Save a figure'''
        plt.savefig(filename, dpi=200)
        
    #ToDo:
    #       *Implement use of missing data coding.
    #       *Provide some way to re-use window (e.g. by
    #         returning w, and having it be an optional
    #         argument. How to clear a window for re-use?
    #       *Make some use of dashed lines in line styles
    #       *The behavior of axis options like reverseX and
    #        XlogAxis is confusing when switchXY = True. We
    #        need to think about the semantics of what we
    #        mean by the "X" axis, and stick to it.  As
    #        currently implemented the "X" axis means the
    #        horizontal axis for these options.  (However,
    #        the semantics is different for the labeling options!)
    #        If we change this (and we probably should), the
    #        scripts plotting soundings will also need to be changed.
    #        as well as the Workbook intructions and problem sets.

#List of line colors and line styles to use
lineColors = ['b','g','r','c','m','y','k']
lineStyles = ['-', '--', '-.', ':'] 
lineThickness = [2, 3, 4]
plotSymbols = ['.','o','v','<','s','*','+','x']

def plot(*curves):
    '''
    Plots Curve objects. plot(...) takes any number of Curve objects
    as arguments, and plots all of them on a single graph.  If multiple
    Curve objects are passed as arguments, the plot appearance data
    (axis options, title, axis labels, etc) are taken from the first
    Curve object. The Curve objects need not all have the same data length.
    This is convenient for plotting data from multiple sources, with varying
    resolution. 
    '''
    cMaster = curves[0] #Appearance options are taken from this Curve
    
    fig = plt.figure() #Always start a new figure
    #**ToDo: Make sure log-axis options so they work consistently
    #        with Ngl implementation when x/y axes are switched.
    #        (Looks OK, but is that implementation the logical way?)
    #r.trXLog = c.XlogAxis
    #r.trYLog = c.YlogAxis
    if cMaster.XlogAxis:
        plt.semilogx()
    if cMaster.YlogAxis:
        plt.semilogy()
    if cMaster.XlogAxis & cMaster.YlogAxis:
        plt.loglog()
    
    plt.title(cMaster.PlotTitle)
    
    #Axis labels (ToDo: add defaults)
    if cMaster.switchXY:
        plt.ylabel(cMaster.Xlabel)
        plt.xlabel(cMaster.Ylabel)
    else:
        plt.xlabel(cMaster.Xlabel)
        plt.ylabel(cMaster.Ylabel)
        
    #  Legends, for multicurve plot
    legends = []
    for c in curves:
        for id in c.listVariables():
            if not id == c.Xid:
                if len(c.label[id]) > 0:
                    legends.append(c.label[id])
                else:
                    legends.append(id)
    #ToDo: Add option to skip legends
    #
    #Suppress line drawing and just plot symbol for scatter plot curves
    #**ToDo: Implement for MatPlotLib
    #r.xyMarkers = plotSymbols
    #r.xyMarkerColors = lineColors
    #r.xyMarkerSizeF = .01
    #r.xyMarkLineModes = []
    formatList = []
    count = 0
    for c in curves:
        for id in c.listVariables():
            if not id == c.Xid:
                if c.scatter[id]:
                    #r.xyMarkLineModes.append('Markers')
                    color = lineColors[count%len(lineColors)]
                    symbol = plotSymbols[count%len(plotSymbols)]
                    formatList.append(color+symbol) 
                else:
                    color = lineColors[count%len(lineColors)]
                    style = lineStyles[count%len(lineStyles)]
                    formatList.append(color+style)
                count += 1
    
    plotList = []
    #Mainly so we can add legends. Could be done with label = ...
    #Note that plt.plot returns a list of lines, even if there is only
    #     one line in the plot, so for legends to work correctly,
    #     we need to extract just the first element.
    countAll = 0
    for c in curves:
        if cMaster.switchXY:
            count = 0
            for data in c.Y():
                plotList.append(plt.plot(data,c.X(),formatList[countAll])[0])
                count += 1
                countAll += 1
        else:
            count = 0
            for data in c.Y():
                plotList.append(plt.plot(c.X(),data,formatList[countAll])[0])
                count += 1
                countAll += 1
                
    #Do the legends
    #11/12/2013: Added Jonah's trick to put legends on the side
    ax = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    plt.legend(plotList, legends, loc='center left', bbox_to_anchor=(1, 0.5))
    
    #Do the axis reversal
    #We do it here, since we don't know the axis limits until plotting is done
    #r.trYReverse = c.reverseY
    #r.trXReverse = c.reverseX
    axes = plt.gca() #Gets current axis
    if cMaster.reverseX:
        axes.set_xlim(axes.get_xlim()[::-1]) #::-1 reverses the array
    if cMaster.reverseY:
        axes.set_ylim(axes.get_ylim()[::-1])
        
    #Now re-draw the plot
    plt.draw()
    
    #(Insert commands needed to show plot, if necessary)
    
    #Eventually we will use this to make subplots and do save option
    return plotObj(None,fig)

def contour(A,**kwargs):
    '''
    A basic contour plotter, which will plot a contour plot
    of a Numeric array. The x and y scales can optionally
    be specified using keyword arguments x and y. For example,
    if we want the x scale to be the array (or list) lat, and
    the y scale to be the array (or list) lon, we would call
    contour as contour(A,x=lat,y=lon).
    '''
    
    #**ToDo: Add labeled contour lines, option to change contour levels,
    #        and to change palette.
    
    fig = plt.figure() #Always start a new figure
    
    # Set axes if they have been specified as keyword arguments
    if 'x' in list(kwargs.keys()):
        x = kwargs['x']
    else:
        x = list(range(A.shape[1]))
    if 'y' in list(kwargs.keys()):
        y = kwargs['y']
    else:
        y = list(range(A.shape[0]))
    cs = plt.contourf(x, y, A)
    cbar = pl.colorbar(cs)
    
    return plotObj(None,fig)
