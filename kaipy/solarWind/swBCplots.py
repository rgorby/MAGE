"""
Generate plots of 1d time series data stored as pyLTR.TimeSeries
objects.  See examples/TimeSeriesPlots.py for usage.
"""
import pylab as p
from kaipy.solarWind.TimeSeries import TimeSeries
import datetime
from matplotlib import dates

def BasicPlot(VarDict,Xname,Yname,Xlabel=True,color='b'):
    """
    Mostly a wrapper for pylab.plot(...)
    """
    x=VarDict[Xname]
    y=VarDict[Yname]
        
    # y may be a time series of tuples/lists of variables to plot 
    # simultaneously (e.g., vector components) using different colors
    if (p.array(y['data'][0]).size > 1 and 
        p.array(y['data'][0]).size == len(color) and
        all([ylen == p.array(y['data'][0]).size for ylen in map(len,y['data'])]) ):
       for i in range(p.array(y['data'][0]).size):
          p.plot(x['data'], [yts[i] for yts in y['data']], color=color[i])
    else:
       p.plot(x['data'],y['data'],color=color)
    
    # Xname may point to a list of datetime.datetime objects, in which case
    # pyplot must be told to plot these as datetimes
    if all([type(ts) == datetime.datetime for ts in x['data']]):
       dfmt = dates.DateFormatter('%m/%d/%y-%H:%M')
       p.gca().xaxis.set_major_formatter(dfmt)
    elif any([type(ts) == datetime.datetime for ts in x['data']]):
       raise Exception('Cannot mix datetime time-axis elements with other types')
    
    
    if Xlabel:
        #locs,labels=p.xticks()
        xStr = x['name']
        if len(x['units']) > 0:
            xStr += ' ['+x['units']+']'
        p.xlabel(xStr)
    else: 
        p.xlabel(' ')
        #locs,labels=p.xticks()
        #p.xticks(locs,(' '))
    
    # y['name'] may not be a scalar
    if (p.array(y['name']).size) > 1:
      p.ylabel('['+y['units']+']',fontsize='small')
    else:
      p.ylabel(y['name']+' ['+y['units']+']',fontsize='small')
  
def SummaryPlot(VarDict,Xname):
    """
    Plot every variable in VarDict.

    This is a simple wrapper around the more generic MultiPlotN.  This
    code may have problems with varDicts that store non-time series
    data.  You've beeen warned.
    """
    #Loop over elements in Dict to make sure they all have data 
    plotVariables=[]
    for var in VarDict.keys():
        if isinstance(VarDict[var], dict):
            if 'data' in VarDict[var]:
               plotVariables.append(var)

    MultiPlotN([VarDict], Xname, plotVariables)


def MultiPlot(VarDict,Xname,Items,color='b'):
    """
    Plot variables stored in TimeSeries object 'varDict'.

    This is a simple wrapper around the more generic MultiPlotN.
    """
    MultiPlotN([VarDict], Xname, Items, [color],[])
    
def MultiPlot2(VarDict,VarDict2,Xname,Items,color1='b',color2='r'):
    """
    Plot items (stored in TimeSeries objects VarDict and VarDict2)
    against one another.  One sub plot for each item.

    This is a simple wrapper around the more generic MultiPlotN.
    """
    MultiPlotN([VarDict, VarDict2],
               Xname,
               Items,
               [color1, color2],
               [])

def MultiPlotN(varDicts, Xname, variables, colors = [], legendLabels=[]):
    """
    Creates one subplot for each variable.  Each subplot renders a plot
    of that variable for each varDict passed.
    
    For example:
      - one varDict and 5 variables will give 5 subplots with one line each
      - 5 varDicts and one variable will give 1 subplot with 5 lines
    
    Parameters:
      varDicts:  List of TimeSeries data dictionaries
      Xname:  key of horizontal axes label (must be defined in all varDicts)
      variables:  List of keys to plot for each varDict
      colors:  color of each line to be drawn
      legendLabels: Display a legend on the plot    
    """
    nSubplots = len(variables)

    # Set default colors (all blue)
    if not colors:
        for d in varDicts:
            colors.append('b')    
  
    # Need to declare these here for proper scope:
    axes=None
    ax=None
  
    for plotIndex, variable in enumerate(variables):
        # Sharing the x-axes applies the same zoom to all subplots when
        # user interacts with a single subplot.
        if plotIndex == 0:
            axes = p.subplot(nSubplots, 1, plotIndex+1)
            ax=axes
        else:
            ax =  p.subplot(nSubplots, 1, plotIndex+1, sharex=axes)
    
        # Turn off x-axis to prevent pointillism problem with upper subplots.
        #ax.xaxis.set_visible(False)
        
        # Alternate vertical axes
        if plotIndex%2 == 0:
            #ax.yaxis.tick_left()
            ax.yaxis.set_ticks_position('left')
            ax.yaxis.set_label_position('left')
        else:
            #ax.yaxis.tick_right()
            ax.yaxis.set_ticks_position('right')
            ax.yaxis.set_label_position('right')
    
        # Fill in subplot data
        for (idx, data) in enumerate(varDicts):
            if plotIndex < nSubplots-1:
                BasicPlot(data, Xname, variable, Xlabel=False, color=colors[idx])
                
                # remove xticklabels from all but bottom subplot...it seems like
                # someone was attempting something similar with ax.xaxis.set_visible(False)
                # in the past, then commented this out; might be worth asking why -EJR 2/2014
                p.setp(ax.get_xticklabels(), visible=False)
            else:
                BasicPlot(data, Xname, variable, Xlabel=True, color=colors[idx])
                
    p.subplots_adjust(hspace=0)
    #ax.xaxis.set_visible(True)
    
    p.subplot(nSubplots, 1, 1)
    if legendLabels:
        p.legend(legendLabels, loc='best')
