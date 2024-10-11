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

from PIDPerfScripts.Definitions import *
from PIDPerfScripts.Binning import *
from PIDPerfScripts.DataFuncs import *
from PIDPerfScripts.PerfCalcFuncs import *
from PIDPerfScripts.RunDictFuncs import *
from PIDPerfScripts.PlotFuncs import *

import ROOT
from ROOT import RooFit


def CalcAveragePerfPlotsAndWriteToFile(TrackDict,
                                       StripVer,
                                       MagPolarity,
                                       PIDCutString,
                                       xBin,
                                       yBin,
                                       zBin,
                                       runMin,
                                       runMax,
                                       schemeName=None,
                                       outputDir=None,
                                       verbose=True,
                                       allowMissingDataSets=False,
                                       maxFiles=-1):

    for k, d in list(TrackDict.items()):
        #print "in calcaverage"
        Plots = GetPerfPlotList(MakePerfPlotsListPyth, StripVer, MagPolarity,
                                d['TrackType'], d['DLLCuts'], d['TrackCuts'],
                                PIDCutString, xBin, yBin, zBin, d['Binning'],
                                runMin, runMax, verbose, allowMissingDataSets,
                                maxFiles)
        #======================================================================
        # Make Weighted Average
        #======================================================================
        MakeAveragePlotList(Plots)
        print('Plots:', Plots)

        #======================================================================
        # Write TH1Fs to file
        #======================================================================
        fileSuffix = ""
        if schemeName is not None:
            fileSuffix += '_{0}'.format(schemeName)
        fileSuffix += '_{0}'.format((d['Binning'][0]).GetName())
        fname = "PubPlots_{plotType}_Strip{strp}_{pol}{suf}.root".format(
            plotType=k, strp=StripVer, pol=MagPolarity, suf=fileSuffix)
        if outputDir is not None:
            fname = "{0}/{1}".format(outputDir, fname)

        f_Out = ROOT.TFile.Open(fname, 'UPDATE')
        if not f_Out:
            raise IOError("Failed to open file {0} for writing".format(fname))

        for subPlots in Plots:
            for plot in subPlots:
                if len(plot) == 1:
                    plot.Write()
                else:
                    for uberplot in plot:
                        uberplot.Write()


# plot.Write()

#for subPlots in Plots:
#    print 'subPlots:' subPlots
#for i in range(0, len(Plots)):
#    for plot in Plots[i]:
#        plot.SetDirectory(f_Out)
#        plot.Write()

        f_Out.Close()

        #======================================================================
        # Delete plots
        #======================================================================
        del Plots


def CalcAveragePerfPlotsAndReturn(TrackDict,
                                  StripVer,
                                  MagPolarity,
                                  runMin,
                                  runMax,
                                  verbose=True,
                                  allowMissingDataSets=False,
                                  maxFiles=-1):

    AveragedPlots = {}

    for k, d in list(TrackDict.items()):

        Plots = GetPerfPlotList(MakePerfPlotsListPyth, StripVer, MagPolarity,
                                d['TrackType'], d['DLLCuts'], d['TrackCuts'],
                                d['Binning'], runMin, runMax, verbose,
                                allowMissingDataSets, maxFiles)

        #======================================================================
        # Make Weighted Average
        #======================================================================
        MakeAveragePlotList(Plots)
        print('Plots:', Plots)

        AveragedPlots[k] = Plots

    return AveragedPlots


def PlotVarsAndReturn(PartName,
                      PlotVars,
                      StripVer,
                      MagPolarity,
                      runMin,
                      runMax,
                      PIDCutString,
                      xBin,
                      yBin,
                      zBin,
                      TrackCuts="",
                      verbose=True,
                      allowMissingDataSets=False,
                      schemeName=None,
                      maxFiles=-1):

    CheckStripVer(StripVer)
    CheckMagPol(MagPolarity)
    CheckPartType(PartName)

    #======================================================================
    # Append runNumber limits to TrackCuts
    #======================================================================
    if None not in (runMin, runMax):
        if TrackCuts != '':
            TrackCuts += ' && '
        TrackCuts += 'runNumber>=' + runMin + ' && runNumber<=' + runMax
    if verbose:
        print('Track Cuts: ', TrackCuts)

    minEntries = 0
    SummedPlots = {}

    for PlotVar in PlotVars:
        SummedPlots[PlotVar] = []

    #======================================================================
    # Loop over all calibration files
    #======================================================================
    #files = GetFiles(StripVer,MagPolarity,PartName,runMin,runMax,maxFiles,verbose)

    #Run 1 file access
    if "Turbo" not in StripVer and "Electron" not in StripVer:
        files = GetFiles(StripVer, MagPolarity, PartName, runMin, runMax,
                         maxFiles, verbose)
    #Run 2 file access (not for 2015-2018 electrons)
    elif "Turbo" in StripVer:
        files = GetWGPFiles(StripVer, MagPolarity, PartName, verbose)
    elif "Electron" in StripVer:
        files = GetElectronFiles(StripVer, MagPolarity, verbose)

    index = 0
    for file in files:
        index += 1
        #DataSet = GetDataSet(StripVer, MagPolarity, PartName, TrackCuts, file, verbose,
        #                    allowMissingDataSets, minEntries)
        DataSet = GetDataSet(StripVer, MagPolarity, PartName, TrackCuts,
                             PIDCutString, xBin, yBin, zBin, file, verbose,
                             allowMissingDataSets, minEntries)
        for PlotVar in PlotVars:
            BinSchema = GetBinScheme(PartName, PlotVar, schemeName)
            if DataSet is not None:
                histname = "%s_%s_%d" % (PartName, PlotVar, index)

                RawHist = ROOT.TH1D(histname, histname, BinSchema.numBins(),
                                    BinSchema.array())
                VarList = ROOT.RooArgList(
                    DataSet.Get_Param(BinSchema.GetName()))
                DataSet.fillHistogram(RawHist, VarList)

                SummedPlots[PlotVar].append(RawHist.Clone(histname))
        if DataSet: DataSet.Delete()

    #======================================================================
    # Sum plots
    #======================================================================

    print('Plots:', SummedPlots)

    for PlotVar in PlotVars:
        MakeSumPlot(SummedPlots[PlotVar])

    print('SummedPlots:', SummedPlots)

    return SummedPlots


def PlotVarAndWriteToFile(PartName,
                          PlotVar,
                          StripVer,
                          MagPolarity,
                          PIDCutString,
                          runMin=None,
                          runMax=None,
                          TrackCuts="",
                          outputDir=None,
                          verbose=True,
                          allowMissingDataSets=False,
                          schemeName=None,
                          maxFiles=-1):
    xBin = ""
    yBin = ""
    zBin = ""
    Plots = PlotVarsAndReturn(PartName, [PlotVar], StripVer, MagPolarity,
                              runMin, runMax, PIDCutString, xBin, yBin, zBin,
                              TrackCuts, verbose, allowMissingDataSets,
                              schemeName, maxFiles)

    #======================================================================
    # Construct file name
    #======================================================================

    fileSuffix = ""
    if schemeName is not None:
        fileSuffix += '_{0}'.format(schemeName)
    #fileSuffix+='_{0}'.format(BinSchema.GetName())
    #fname="PlotCalibDistribution_{particle}_{var}_Strip{strp}_{pol}{suf}.root".format(
    #    particle=PartName, var=PlotVar, strp=StripVer, pol=MagPolarity, suf=fileSuffix)
    fname = "PlotCalibDistribution_{particle}_{var}_Strip{strp}_{pol}.root".format(
        particle=PartName, var=PlotVar, strp=StripVer, pol=MagPolarity)
    if outputDir is not None:
        fname = "{0}/{1}".format(outputDir, fname)

    f_Out = ROOT.TFile.Open(fname, 'UPDATE')
    if not f_Out:
        raise IOError("Failed to open file {0} for writing".format(fname))

    for plot in Plots[PlotVar]:
        plot.Write()

    f_Out.Close()

    del Plots


def GetPerfPlotList(PerfFunc,
                    StripVer,
                    MagPolarity,
                    PartName,
                    DLLCutList,
                    TrackCuts,
                    PIDCutString,
                    xBin,
                    yBin,
                    zBin,
                    BinningScheme=None,
                    runMin=None,
                    runMax=None,
                    verbose=True,
                    allowMissingDataSets=False,
                    maxFiles=-1):

    #print " IN getperfplotlist"
    #**********************************************************************
    CheckStripVer(StripVer)
    CheckMagPol(MagPolarity)
    CheckPartType(PartName)
    print(BinningScheme.size())
    print(BinningScheme.at(0))
    #======================================================================
    # Append runNumber limits to TrackCuts
    #======================================================================
    #if None not in (runMin, runMax):
    #    if TrackCuts!='':
    #        TrackCuts+=' && '
    #    TrackCuts+='runNumber>='+runMin+' && runNumber<='+runMax
    #if verbose:
    #    print 'Track Cuts: ', TrackCuts

    #======================================================================
    # Declare default list of PID plots
    #======================================================================
    Plots = []
    minEntries = 0

    #======================================================================
    # Loop over all calibration subsamples
    #======================================================================

    #Use GetFiles function for Run I data
    if "Turbo" not in StripVer and "Electron" not in StripVer:
        files = GetFiles(StripVer, MagPolarity, PartName, runMin, runMax,
                         maxFiles, verbose)
    elif "Turbo" in StripVer:
        files = GetWGPFiles(StripVer, MagPolarity, PartName, verbose)
    elif "Electron" in StripVer:
        files = GetElectronFiles(StripVer, MagPolarity, verbose)

    for file in files:
        #Now pass in all of the variables used, in order to filter them from the WGP nTuple
        DataSet = GetDataSet(StripVer, MagPolarity, PartName, TrackCuts,
                             PIDCutString, xBin, yBin, zBin, file, verbose,
                             allowMissingDataSets, minEntries)
        if DataSet is not None:
            #======================================================================
            # Run Specific implementation of PerfCalculator
            #======================================================================

            #print "printing the size again"
            #print BinningScheme.size()
            Plots = PerfFunc(
                PartName,
                DataSet,
                DLLCutList,
                BinningScheme,
                Plots,
                verbose=verbose)
            #======================================================================
            # Delete DataSet
            #======================================================================
            #            DataSet.Delete()
            del DataSet  # avoids a weird segmentation violation
    #======================================================================
    # Return list of plots
    #======================================================================
    return Plots


def GetPerfResultList(PerfFunc,
                      StripVer,
                      MagPolarity,
                      PartName,
                      DLLCutList,
                      TrackCuts,
                      runMin=None,
                      runMax=None,
                      verbose=True,
                      allowMissingDataSets=False,
                      maxFiles=-1):
    #**********************************************************************
    CheckStripVer(StripVer)
    CheckMagPol(MagPolarity)
    CheckPartType(PartName)

    #======================================================================
    # Append runNumber limits to TrackCuts
    #======================================================================
    #if None not in (runMin, runMax):
    #    if TrackCuts!='':
    #        TrackCuts+=' && '
    #    TrackCuts+='runNumber>='+runMin+' && runNumber<='+runMax
    #if verbose:
    #    print 'Track Cuts: ', TrackCuts

    #======================================================================
    # Declare default list of PID results
    #======================================================================
    PIDResults = []
    minEntries = 0

    #======================================================================
    # Loop over all calibration subsamples
    #======================================================================
    #  for DataSet in GetDataSets(StripVer, MagPolarity, PartName, TrackCuts,
    #                             runMin, runMax, verbose, allowMissingDataSets):
    files = GetFiles(StripVer, MagPolarity, PartName, runMin, runMax, maxFiles,
                     verbose)
    for file in files:
        DataSet = GetDataSet(StripVer, MagPolarity, PartName, TrackCuts, file,
                             verbose, allowMissingDataSets, 1000)
        if DataSet is not None:

            #======================================================================
            # Run Specific implementation of PerfCalculator
            #======================================================================
            PerfFunc(DataSet, DLLCutList, PIDResults, verbose)
            if verbose:
                print(PIDResults)

        #======================================================================
        # Delete DataSet
        #======================================================================
        del DataSet
        #DataSet.Delete()

    #======================================================================
    # Return list PIDResults
    #======================================================================
    return PIDResults
