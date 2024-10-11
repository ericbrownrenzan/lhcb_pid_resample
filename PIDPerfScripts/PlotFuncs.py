###############################################################################
# (c) Copyright 2021 CERN for the benefit of the LHCb Collaboration           #
#                                                                             #
# This software is distributed under the terms of the GNU General Public      #
# Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   #
#                                                                             #
# In applying this licence, CERN does not waive the privileges and immunities #
# granted to it by virtue of its status as an Intergovernmental Organization  #
# or submit itself to any jurisdiction.                                       #
###############################################################################



from past.utils import old_div
import ROOT
from array import array


def MakeAveragePlot(Plots, oldStyleAverage=False):
    #======================================================================
    # Loop over plots in subPlots and append average plot
    # oldStyleAverage enables a bad way of calculating the average which
    # does not work for small efficiencies, enable only if you want to
    # compare with old PIDCalib versions
    #======================================================================
    CombPlot = None
    SumTotalsPlot = None
    SumPassedPlot = None

    if len(Plots) > 0 and isinstance(Plots[0], list):
        if oldStyleAverage:
            print("")
            print(
                "Warning: using old style of average calculation. Efficiencies determined this way are incorrect if the efficiency is near zero."
            )
            print("")
            for (average, total, passed) in Plots:
                average.SetBit(ROOT.TH1.kIsAverage)
                if CombPlot is None:
                    CombPlot = average.Clone(average.GetName().replace(
                        '_0', '_All'))
                    SumTotalsPlot = total.Clone(total.GetName().replace(
                        '_0', '_All'))
                    SumPassedPlot = passed.Clone(passed.GetName().replace(
                        '_0', '_All'))
                else:
                    CombPlot.Add(average)
                    SumTotalsPlot.Add(total)
                    SumPassedPlot.Add(passed)
            Plots.append([CombPlot, SumTotalsPlot, SumPassedPlot])
        else:
            print(" Using new style average ")
            for (average, total, passed) in Plots:
                #plot.SetBit(ROOT.TH1.kIsAverage)
                if CombPlot is None:
                    #print plot.GetName()
                    CombPlot = average.Clone(average.GetName().replace(
                        '_0', '_All'))
                    SumTotalsPlot = total.Clone(total.GetName().replace(
                        '_0', '_All'))
                    SumPassedPlot = passed.Clone(passed.GetName().replace(
                        '_0', '_All'))
                else:
                    SumTotalsPlot.Add(total)
                    SumPassedPlot.Add(passed)
            CombPlot.Divide(SumPassedPlot, SumTotalsPlot, 1.0, 1.0, "B")
            Plots.append([CombPlot, SumTotalsPlot, SumPassedPlot])
    else:
        for plot in Plots:
            plot.SetBit(ROOT.TH1.kIsAverage)
            if CombPlot is None:
                #print plot.GetName()
                CombPlot = plot.Clone(plot.GetName().replace('_0', '_All'))
            else:
                CombPlot.Add(plot)
        Plots.append(CombPlot)
    #print Plots


def MakeSumPlot(Plots):
    #======================================================================
    # Loop over plots in subPlots and append sum plot
    #======================================================================
    CombPlot = None
    for plot in Plots:
        if CombPlot is None:
            #print plot.GetName()
            CombPlot = plot.Clone(plot.GetName().replace('_0', '_All'))
        else:
            CombPlot.Add(plot)
    Plots.append(CombPlot)
    #print Plots


def MakeAveragePlotList(Plots, oldStyleAverage=False):
    #======================================================================
    # Loop over all subPlots and run MakeAveragePlot()
    #======================================================================
    for i in range(0, len(Plots)):
        MakeAveragePlot(Plots[i], oldStyleAverage)


def MakeSumPlotList(Plots):
    #======================================================================
    # Loop over all subPlots and run MakeSumPlot()
    #======================================================================
    for i in range(0, len(Plots)):
        MakeSumPlot(Plots[i])


def CalculateBinning(plot, nBins):
    #tlines = []
    points = []
    totalIntegral = plot.Integral()
    plotBins = plot.GetNbinsX()
    for i in range(nBins - 1):
        j = 1
        integral = 0.
        xPos = None
        for j in range(1, plotBins + 2):
            integral += plot.GetBinContent(j)
            #print i,j,integral,(i,1.),totalIntegral,nBins
            if integral > old_div((i + 1.) * totalIntegral, nBins):
                xPos = plot.GetBinCenter(j)
                break
            j += 1
        points.append(xPos)
        #tlines.append(TLine(xPos,1.1,xPos,1.))
        #tlines[-1].SetLineColor(colour)
    return points


def makeTGraph(PIDResults, MisIDResults, name, title=None):
    #======================================================================
    # Create TGraph
    #======================================================================
    x, x_err = [t for t in zip(*PIDResults)]
    y, y_err = [t for t in zip(*MisIDResults)]
    ## x     = [list(t) for t in zip(*KaonResults)][0]
    ## x_err = [list(t) for t in zip(*KaonResults)][1]
    ## y     = [list(t) for t in zip(*PionResults)][0]
    ## y_err = [list(t) for t in zip(*PionResults)][1]

    graph = ROOT.TGraphErrors(
        len(x), array('f', x), array('f', y), array('f', x_err),
        array('f', y_err))
    graph.SetName(name)
    if title is not None:
        graph.SetTitle(title)
    return graph
