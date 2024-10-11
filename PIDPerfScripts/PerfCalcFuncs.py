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
from __future__ import print_function
from __future__ import division

from past.utils import old_div
import ROOT
from ROOT import RooFit
from ROOT import vector
from operator import add
from math import sqrt
import warnings
import io
import sys
from PIDPerfScripts.Definitions import *
from PIDPerfScripts.Binning import GetBinScheme


def MakePerfPlotsList(PartName,
                      DataSet,
                      DLLCutList,
                      BinningScheme,
                      Plots=None,
                      verbose=True):

    #======================================================================
    # Declare instances of PerfCalculator tool : (class)(instance of class)
    #======================================================================
    ROOT.gSystem.Load('libRooStats.so')
    #    ROOT.gSystem.Load('libCintex.so')
    #   cintex=ROOT.Cintex
    #  cintex.Enable()
    ROOT.gSystem.Load('libPIDPerfToolsLib.so')
    ROOT.gSystem.Load('libPIDPerfToolsDict.so')
    Calc = ROOT.PerfCalculator(ROOT.GenericDataSet)(DataSet)

    #======================================================================
    # Define array of lists, one for each DLL cut
    #======================================================================
    subPlots = [[] for j in range(len(DLLCutList))]

    #======================================================================
    # Loop over DLL cuts, and append TH1F to appropriate subPlot-list
    #======================================================================
    ROOT.TH1.SetDefaultSumw2(True)
    dataset_index = len(Plots[0]) if Plots else 0
    for index, cut in enumerate(DLLCutList):
        parsed_cut = cut.strip().replace("/", "_div_")
        histname = "%s_%s_%d" % (PartName, parsed_cut, dataset_index)
        Dummy = Calc.Perf_Hist(histname, cut, BinningScheme, False, 5, verbose)
        subPlots[index].append(Dummy)
        #Dummy.Delete()
    if verbose:
        print('subPlots: ', subPlots)
    del Calc

    #======================================================================
    # Append subPlots to Plots, i.e.
    # Plots = [ [0,1,...]_DLLK>A, [0,1,...]_DLLK>B, ... ]
    #======================================================================
    # Need to resolve the following line....
    NewPlots = list(map(add, Plots, subPlots)) if Plots else subPlots
    if verbose:
        print('NewPlots: ', dataset_index, NewPlots)
    del Plots

    return NewPlots


def MakePerfPlotsListPyth(PartName,
                          DataSet,
                          DLLCutList,
                          BinningScheme,
                          Plots=None,
                          verbose=True):
    #print "in perfplots list pyth"
    ROOT.TH1.SetDefaultSumw2(True)
    dataset_index = len(Plots[0]) if Plots else 0
    PlotsByCut = []
    #print " abit further"
    print(BinningScheme.at(0))

    if BinningScheme.size() > 0:
        BinSchema1 = BinningScheme.at(0)
    if BinningScheme.size() > 1:
        BinSchema2 = BinningScheme.at(1)
    if BinningScheme.size() > 2:
        BinSchema3 = BinningScheme.at(2)

    #print "understood bin scheme"
    for index, cut in enumerate(DLLCutList):
        parsed_cut = cut.strip().replace("/", "_div_")
        histname = "%s_%s_%d" % (PartName, parsed_cut, dataset_index)

        #   print "understood the Cut list"
        YVarArgVar = ROOT.RooCmdArg.none()
        ZVarArgVar = ROOT.RooCmdArg.none()
        if BinningScheme.size() > 1:
            YVarArgVar = RooFit.YVar(
                DataSet.Get_Param(BinSchema2.GetName()),
                RooFit.Binning(BinSchema2))
        if BinningScheme.size() > 2:
            ZVarArgVar = RooFit.ZVar(
                DataSet.Get_Param(BinSchema3.GetName()),
                RooFit.Binning(BinSchema3))

        TotalHist = ROOT.RooAbsData.createHistogram(
            DataSet, "TotalHist_" + histname,
            DataSet.Get_Param(BinSchema1.GetName()),
            RooFit.Binning(BinSchema1), YVarArgVar, ZVarArgVar)

        PassedHist = ROOT.RooAbsData.createHistogram(
            DataSet, "PassedHist_" + histname,
            DataSet.Get_Param(BinSchema1.GetName()),
            RooFit.Binning(BinSchema1),
            RooFit.Cut(DataSet.ReplaceVarNamesInCuts(cut)), YVarArgVar,
            ZVarArgVar)

        if BinningScheme.size() > 2:
            RatioHist = ROOT.TH3F(histname, histname, BinSchema1.numBins(),
                                  BinSchema1.array(), BinSchema2.numBins(),
                                  BinSchema2.array(), BinSchema3.numBins(),
                                  BinSchema3.array())
            BinSchemas = [BinSchema1, BinSchema2, BinSchema3]
            AxisTitles = [
                RatioHist.SetXTitle, RatioHist.SetYTitle, RatioHist.SetZTitle
            ]
        elif BinningScheme.size() > 1:
            RatioHist = ROOT.TH2F(histname, histname, BinSchema1.numBins(),
                                  BinSchema1.array(), BinSchema2.numBins(),
                                  BinSchema2.array())
            BinSchemas = [BinSchema1, BinSchema2]
            AxisTitles = [RatioHist.SetXTitle, RatioHist.SetYTitle]
        elif BinningScheme.size() > 0:
            RatioHist = ROOT.TH1F(histname, histname, BinSchema1.numBins(),
                                  BinSchema1.array())
            BinSchemas = [BinSchema1]
            AxisTitles = [RatioHist.SetXTitle]

        for scheme, setTitle in zip(BinSchemas, AxisTitles):
            title = DataSet.Get_Param(scheme.GetName()).GetTitle()
            unit = DataSet.Get_Param(scheme.GetName()).getUnit().strip()
            if unit != "":
                title += " [%s]" % (unit)

            # print "-----"
            # print setTitle
            # DataSet.Get_Param(scheme.GetName()).Print()
            # print title
            # print unit
            # print scheme.GetName()
            setTitle(title)
            # print "-----"

        #RatioHist = PassedHist.Clone(histname)
        RatioHist.Divide(PassedHist, TotalHist, 1.0, 1.0, "B")
        #print "done some dividing"
        if BinningScheme.size() > 2:
            pidTable = ROOT.PIDTable(RatioHist, BinSchema1, BinSchema2,
                                     BinSchema3)
        elif BinningScheme.size() > 1:
            pidTable = ROOT.PIDTable(RatioHist, BinSchema1, BinSchema2)
        elif BinningScheme.size() > 0:
            pidTable = ROOT.PIDTable(RatioHist, BinSchema1)
        #print "pidTable made"
        #pidTable.PrintTable(sys.stdout, False)
        #sys.stdout.flush()
        #print "but not printed the table"
        PlotsByCut.append([[RatioHist, TotalHist, PassedHist]])

    NewPlots = list(map(add, Plots, PlotsByCut)) if Plots else PlotsByCut
    if verbose:
        print('PlotsByCut: ', PlotsByCut)

    return NewPlots


def MakePIDResultList(DataSet, DLLCutList, PIDResults, verbose=True):
    print(PIDResults)
    #======================================================================
    # Declare instances of PerfCalculator tool : (class)(instance of class)
    #======================================================================

    ROOT.gSystem.Load('libRooStats.so')
    #ROOT.gSystem.Load('libCintex.so')
    #cintex=ROOT.Cintex
    #cintex.Enable()
    ROOT.gSystem.Load('libPIDPerfToolsLib.so')
    ROOT.gSystem.Load('libPIDPerfToolsDict.so')

    Calc = ROOT.PerfCalculator(ROOT.GenericDataSet)(DataSet)

    #======================================================================
    # Calculate Efficiency
    #======================================================================
    Dummy = Calc.Perf(DLLCutList, False)  #, verbose)
    PIDResults.append(Dummy)

    print(PIDResults)


def WeightedAveragePIDResults(PIDResults, verbose=False):
    #======================================================================
    # Loop over list of PIDResult::Container (i.e. vector<PIDResult>) and
    # perform weighted average for each index, returning a list of
    # [[Eff_0, Err_0], [Eff_1, Err_1], ... , [Eff_n, Err_n]]
    #======================================================================
    nresults = PIDResults[0].size()
    CombPIDResult = []
    if nresults == 0:
        msg = ("Failed to perform weighted average of PID results. "
               "No results in PIDResult::Container.")
        warnings.warn(msg)
        return CombPIDResult

    for result in PIDResults:
        if result.size() != nresults:
            raise ValueError(
                'PIDResult::Container objects are of differing size')

    for ieff in range(0, nresults):
        #print eff
        Eff_Num = 0.0
        Eff_Denom = 0.0
        ##Err_InSq = 0.0
        for iresult, result in enumerate(PIDResults):
            eff = result[ieff].Eff
            err = result[ieff].Err_Lo
            if err == 0.:
                msg = (
                    "Got PIDResult::Container result (index {0:d}, result {1:d}) "
                    "with zero error. This result will be excluded from the "
                    "weighted average").format(ieff, iresult)
                warnings.warn(msg)
                continue
            Eff_Num += old_div((1.0 * eff), (err * err))
            Eff_Denom += 1.0 / (err * err)
            ##Err_InSq += 1.0/(err*err)
        if Eff_Denom == 0.0:
            msg = (
                "Got PIDResult::Container element (index {0:d}) "
                "with zero error. This index will be excluded from the list of "
                "averaged PID results").format(ieff)
            warnings.warn(msg)
            continue
        CombPIDResult.append(
            [old_div(Eff_Num, Eff_Denom),
             sqrt(1.0 / Eff_Denom)])
        if verbose:
            msg = ("PIDResult::Container index {idx:d}"
                   ", weighted efficiency=({eff:.3f}+/-{err:.3f})%").format(
                       idx=ieff, eff=eff * 100.0, err=err * 100.0)
            print(msg)
            ##sqrt(1/Err_InSq)])
    #print CombPIDResult

    return CombPIDResult


def SimpleAveragePIDResults(CombPIDResults_a, CombPIDResults_b):
    #======================================================================
    # Perform the simple average of CombPIDResults_a[i] with
    # CombPIDResults_b[i] and return another list of type
    # [[Eff_0, Err_0], [Eff_1, Err_1], ... , [Eff_n, Err_n]]
    #======================================================================
    if len(CombPIDResults_a) != len(CombPIDResults_b):
        print('**** Efficiency lists of differing size ****')
        sys.exit(1)

    CombPIDResult = []
    for i in range(0, len(CombPIDResults_a)):
        Eff = old_div((CombPIDResults_a[i][0] + CombPIDResults_b[i][0]), 2)
        Err_Sq = pow(old_div(CombPIDResults_a[i][1], 2), 2) + pow(
            old_div(CombPIDResults_b[i][1], 2), 2)
        CombPIDResult.append([Eff, sqrt(Err_Sq)])

    return CombPIDResult
