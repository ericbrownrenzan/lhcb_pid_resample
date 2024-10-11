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
from PIDPerfScripts.Definitions import *
from PIDPerfScripts import OverrideCalibDataStore

#Make itertool.zip_longest usable for Python 2 and 3 following:
#https://python-future.org/compatible_idioms.html
from future.moves.itertools import filterfalse, zip_longest
import pickle
import os
import numpy as np


def GetPklFileSuffix(PartName):
    #print "Getting pkl file suffix "
    #print PartName
    CheckPartType(PartName)
    #print "passed checkplktype "
    if PartName in ("K", "Pi"):
        return "dst_k_and_pi"
    if PartName in ("P", ):
        return "lam0_p"
    if PartName in ("Mu", ):
        return "jpsi_mu"
    if PartName in ("e", ):
        return "jpsi_e"
    if PartName in ("K_MuonUnBiased", "Pi_MuonUnBiased"):
        return "dst_k_and_pi_muonUnBiased"
    if PartName in ("P_MuonUnBiased", ):
        return "lam0_p_muonUnBiased"
    if PartName in ("e_MuonUnBiased", ):
        return "jpsi_e_muonUnBiased"
    if PartName in ("P_LcfB", ):
        return "LcfB_P"
    if PartName in ("P_TotLc", ):
        return "TotLc_P"
    if PartName in ("P_IncLc", ):
        return "IncLc_P"


def GetRunDictionary(StripVer, PartName="K",
                     verbose=True):  #, IsMuonUnBiased=False):

    pklfileSuffix = GetPklFileSuffix(PartName)
    #MUONPreFix = ''
    ##   if PartName in ("Mu", "K_MuonUnBiased", "Pi_MuonUnBiased", "P_MuonUnBiased"):
    ##         MUONPreFix = '/MUON'

    #======================================================================
    # Dictionary of Dictionaires for StripVersion -> {UpRuns, DownRuns}
    #======================================================================
    stv = StripVer

    if StripVer == '21_MCTuneV4':
        stv = '21'
    if StripVer == '21r1_MCTuneV4':
        stv = '21r1'
    if StripVer == '23_MCTuneV1':
        stv = '23'

    UpRunLims = pickle.load(
        open(
            os.path.expandvars(
                '$PIDPERFSCRIPTSROOT/pklfiles/Stripping{strp}/up_runLimits_{suf}.pkl'
                .format(strp=stv, suf=pklfileSuffix)), 'rb'),
        encoding='latin1')
    DownRunLims = pickle.load(
        open(
            os.path.expandvars(
                '$PIDPERFSCRIPTSROOT/pklfiles/Stripping{strp}/down_runLimits_{suf}.pkl'
                .format(strp=stv, suf=pklfileSuffix)), 'rb'),
        encoding='latin1')

    StripDict = {
        'StripVer': StripVer,
        'RecoVer': GetRecoVer(StripVer),
        'UpRuns': UpRunLims.astype(int),
        'DownRuns': DownRunLims.astype(int)
    }

    NumberOfRunColumns = 7

    UpRunsStr = zip_longest(
        *[iter([l.__repr__()
                for l in StripDict['UpRuns'].tolist()])] * NumberOfRunColumns,
        fillvalue="")
    UpRunsStr = [
        ", ".join([s for s in arr if s != ""]) + "," for arr in UpRunsStr
    ]

    DownRunsStr = zip_longest(
        *[iter([l.__repr__() for l in StripDict['DownRuns'].tolist()])] *
        NumberOfRunColumns,
        fillvalue="")
    DownRunsStr = [
        ", ".join([s for s in arr if s != ""]) + "," for arr in DownRunsStr
    ]

    if verbose:

        print('StripVer   : ', StripDict['StripVer'])
        print('RecoVer    : ', StripDict['RecoVer'])
        print('UpRuns     :  [')
        for s in UpRunsStr:
            print("    " + s)
        print("]")
        print('DownRuns   :  [')
        for s in DownRunsStr:
            print("    " + s)
        print("]")

    return StripDict


def __FindFileIndex(run, runFirst, runLast):
    #======================================================================
    # Is run within the limits of a file?
    # - Create arrays of booleans determining if run is:
    #    a) Equal-to-or-greater-than the first run in each
    #       file (IsrunGreaterThanFirst)
    #    b) Equal-to-or-less-than the last run in each
    #       file (IsrunLessThanLast)
    # - Take the logical and and search for a 'True'.
    # - If there is no 'True', then find the first 'True' in
    #   IsrunLessThanLast.
    # - If there is no 'True' in IsrunLessThanLast, then find the
    #   first 'True' in IsrunGreaterThanFirst.
    #   and take this index as the first file to process.
    # - If no match is made, a 'ValueError: min() arg is an empty sequence'
    #   is returned.
    #======================================================================
    IsrunGreaterThanFirst = runFirst <= run
    IsrunLessThanLast = runLast >= run
    IsrunInHere = np.logical_and(IsrunGreaterThanFirst, IsrunLessThanLast)

    if np.any(IsrunInHere):
        FileIndex = np.min(np.nonzero(IsrunInHere == True)[0])
    elif np.any(IsrunLessThanLast):
        FileIndex = np.min(np.nonzero(IsrunLessThanLast == True)[0])
    else:
        FileIndex = np.max(np.nonzero(IsrunGreaterThanFirst == True)[0])

    return FileIndex


def GetMinMaxFileDictionary(DataDict,
                            MagPolarity,
                            runMin=None,
                            runMax=None,
                            maxFiles=-1,
                            verbose=True):
    overridden = OverrideCalibDataStore.GetDictIndex(runMin, runMax, maxFiles,
                                                     verbose)
    if overridden: return overridden
    #======================================================================
    # Determine file index ranges corresponding to RunMin and RunMax
    #======================================================================
    CheckMagPol(MagPolarity)
    runFirst = 0
    runLast = 0
    if MagPolarity == 'MagUp':
        ndim = DataDict['UpRuns'].shape[0]
        runFirst = DataDict['UpRuns'][0:ndim, 0]
        runLast = DataDict['UpRuns'][0:ndim, 1]
    else:
        ndim = DataDict['DownRuns'].shape[0]
        runFirst = DataDict['DownRuns'][0:ndim, 0]
        runLast = DataDict['DownRuns'][0:ndim, 1]

    if runMin is None:
        #print 'runMin is None', runMin
        minIndex = 0
    else:
        #print 'runMin is not None', runMin
        minIndex = __FindFileIndex(int(runMin), runFirst, runLast)
    if runMax is None:
        #print 'runMax is None', runMax
        maxIndex = len(runFirst) - 1
    else:
        #print 'runMax is not None', runMax
        maxIndex = __FindFileIndex(int(runMax), runFirst, runLast)

    #print 'runFirst: {runFirst}, runLast: {runLast}'.format(runFirst=runFirst, runLast=runLast)
    if verbose:
        print('RunMin: {run}, Min file index: {index}'.format(
            run=runMin, index=minIndex))
        print('RunMax: {run}, Max file index: {index}'.format(
            run=runMax, index=maxIndex))

    if maxFiles != -1 and maxFiles != None:
        if int(maxFiles) - 1 < (maxIndex - minIndex):
            maxIndex = minIndex + int(maxFiles) - 1
            if verbose:
                print('MaxFiles: {run}, Min file index: {index}'.format(
                    run=maxFiles, index=minIndex))
                print('MaxFiles: {run}, Max file index: {index}'.format(
                    run=maxFiles, index=maxIndex))

    ret = {'minIndex': minIndex, 'maxIndex': maxIndex}

    return ret
