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
import ROOT
import sys
import warnings
import time
import os
from os import path
from PIDPerfScripts.Definitions import *
from PIDPerfScripts.Exceptions import *
from PIDPerfScripts.RunDictFuncs import *
from PIDPerfScripts.TupleDataset import *
from PIDPerfScripts import OverrideCalibDataStore
import sys


def GetDataSetNameDictionary(PartName):
    #======================================================================
    # Define Mother and Workspace name given Particle name
    #======================================================================
    MotherName = GetMotherName(PartName)
    WSName = GetWorkspaceName(PartName)

    ret = {
        'MotherName': GetMotherName(PartName),
        'WorkspaceName': GetWorkspaceName(PartName)
    }

    return ret


def GetDataSets(StripVer,
                MagPolarity,
                PartName,
                TrackCuts,
                runMin=None,
                runMax=None,
                verbose=False,
                allowMissingDataSets=False,
                minEntries=0,
                maxFiles=-1):

    CheckStripVer(StripVer)
    if None not in (runMin, runMax):
        if TrackCuts != '':
            TrackCuts += ' && '
        TrackCuts += 'runNumber>=' + runMin + ' && runNumber<=' + runMax
    if verbose:
        print('Track Cuts: ', TrackCuts)

    if verbose:
        if sys.version_info[0] < 3:
            print("Attemptng to get URLS: ({0},{1}) ".format(
                time.time(), time.clock()))
        else:
            print("Attemptng to get URLS: ({0},{1}) ".format(
                time.time(), time.process_time()))
    files = GetFiles(StripVer, MagPolarity, PartName, runMin, runMax, maxFiles,
                     verbose)

    if verbose:
        if sys.version_info[0] < 3:
            print("Obtain URLS: ({0},{1}) ".format(time.time(), time.clock()))
        else:
            print("Obtain URLS: ({0},{1}) ".format(time.time(),
                                                   time.process_time()))
    DataSets = []

    for file in files:
        DataSet = GetDataSet(StripVer, MagPolarity, PartName, TrackCuts, file,
                             verbose, allowMissingDataSets, minEntries)
        if DataSet is not None:
            DataSets.append(DataSet)
    return DataSets


def GetDataSet(StripVer,
               MagPolarity,
               PartName,
               TrackCuts,
               PIDCutString,
               xBin,
               yBin,
               zBin,
               file,
               verbose=False,
               allowMissingDataSets=False,
               minEntries=1):

    #If Run I data has been requested and MC12TuneV4 or MC15TuneV1 ProbNN cut used, return error
    if 'Turbo' not in StripVer and 'Electron' not in StripVer:
        if 'MC12TuneV4' in PIDCutString or 'MC15TuneV1' in PIDCutString:
            raise RooWorkspaceError(
                "Cannot use MC12TuneV4 or MC15TuneV1 ProbNN for Run 1 data. Please use MC12TuneV2 or MC12TuneV3"
            )

    RecoVer = GetRecoVer(StripVer)

    if verbose:
        if sys.version_info[0] < 3:
            print("Attempting to open file {0} for reading: ({1},{2})".format(
                file, time.time(), time.clock()))
        else:
            print("Attempting to open file {0} for reading: ({1},{2})".format(
                file, time.time(), time.process_time()))
    f = ROOT.TFile.Open(file)

    if not f:
        if allowMissingDataSets:
            warnings.warn(
                "File %s does not exist. Skipping this subsample" % file)
            return None
        else:
            raise IOError("Failed to open file %s for reading" % file)

    DataSetNameDict = GetDataSetNameDictionary(PartName)
    wsname = DataSetNameDict['WorkspaceName']
    if verbose:
        print("Attempting to get workspace {0}".format(wsname))

    try:
        ws = f.Get(wsname)
    except:
        print("Converting Run 2 data to RooDataSet...")
        ws = None

    #In the case of Run II data, read the WGP nTuple and convert to RooDataSet
    CheckPartType(PartName)
    CheckStripVerPartNameMagPol(StripVer, PartName, MagPolarity)
    PartType = GetRealPartType(PartName)
    DataSetNameDict = GetDataSetNameDictionary(PartName)

    if verbose:
        if sys.version_info[0] < 3:
            print("Attempting to create dataset: ({0},{1})".format(
                time.time(), time.clock()))
        else:
            print("Attempting to create dataset: ({0},{1})".format(
                time.time(), time.process_time()))
    if ws:
        Data = ws.data('data')
    else:
        Data = getDataSetFromTuple(
            file=file,
            mother=DataSetNameDict['MotherName'],
            part=PartType,
            trackcuts=TrackCuts,
            pidcuts=PIDCutString,
            xvar=xBin,
            yvar=yBin,
            zvar=zBin,
            strip=
            StripVer  #Added to check the nTuple names for a corresponding WGP sample
        )
        print("Automatic conversion from TTree to RooDataSet")
    if Data is None:
        print("No entires in RooDataSet - skipping this input file")
        return None

#        raise TFileError("Failed to retrieve workspace {wsname} from file {fname}".format(
#            wsname=DataSetNameDict['WorkspaceName'], fname=fname))
    if not Data:
        raise RooWorkspaceError(
            "RooDataSet not found in workspace %s" % wsname)
    if verbose:
        if sys.version_info[0] < 3:
            print("Retrieved data: ({0},{1})".format(time.time(),
                                                     time.clock()))
        else:
            print("Retrieved data: ({0},{1})".format(time.time(),
                                                     time.process_time()))
    Data.Print("v")

    if verbose:
        print("Applying any cuts to the dataset.")

    #======================================================================
    # Declare Instance of RICHTrackDataSet for Calibration tracks
    #======================================================================
    ROOT.gSystem.Load('libRooStats.so')
    #ROOT.gSystem.Load('libCintex.so')
    #cintex=ROOT.Cintex
    #cintex.Enable()
    ROOT.gSystem.Load('libPIDPerfToolsLib.so')
    ROOT.gSystem.Load('libPIDPerfToolsDict.so')

    DataSet = None
    dsType = None
    if (RecoVer >= 14):
        dsType = 'GenericDataSet'

        VariableVector = ROOT.std.vector(ROOT.std.pair("string,string"))()
        for VarName, DataSetVarName in DataSetVariables().items():
            #VariableAlias = ROOT.std.pair("string,string")(VarName, DataSetVarName.format(particle=PartType))
            VariableAlias = ROOT.std.pair("string,string")(
                VarName, DataSetVarName.dsname.format(particle=PartType))
            VariableVector.push_back(VariableAlias)

        DataSet = ROOT.GenericDataSet('Calib_Data', Data, Data.get(),
                                      VariableVector, TrackCuts, 'nsig_sw')

    else:
        print(
            'Reco Version < 14! Please use an earlier version of PIDPerfTools.'
        )

    if ws:
        ws.Delete()

    f.Close()
    if verbose:
        print("DataSet is", str(DataSet) + ", which uses the internal store:",
              str(DataSet.store()))

    #======================================================================
    # Sanity test: do we have a dataset, and is it empty?
    #======================================================================
    if DataSet is None:
        raise RooDataSetError(
            "Failed to create {0} from RooDataSet".format(dsType))

    if DataSet.sumEntries() == 0:
        raise RooDataSetError("{0} contains no entries".format(dsType))

    if verbose:
        DataSet.Print('v')

    #======================================================================
    # Reduce dataset to only those events within binning limits
    #======================================================================
    #DataSet = None
    #if BinningScheme is not None:
    #    DataSet = AllDataSet.SetInBinSchema(BinningScheme)
    #    AllDataSet.Delete()
    #else:
    #    DataSet = AllDataSet

    #======================================================================
    # Veto ranges with insufficient statistics
    #======================================================================
    if DataSet.sumEntries() < minEntries:
        msg = ('Insufficient events in {dsType} ({nEvt}). '
               'Requested at least {minEvts} entries. '
               'Skipping this sumsample').format(
                   dsType=dsType,
                   nEvt=DataSet.sumEntries(),
                   minEvts=minEntries)
        warnings.warn(msg)
        return None

    return DataSet


##############################################
### RUN II FILE ACCESS WITH THIS FUNCTION ####
##############################################


def GetWGPFiles(StripVer, MagPolarity, PartName, verbose=False):

    import os
    import subprocess
    from subprocess import Popen, PIPE
    import tempfile
    import json

    #Run the Dirac script using the required compiler version for Dirac
    #Pass it a temporary directory to write the file list to
    path = tempfile.mkdtemp()

    arguments = [
        'unset PYTHONHOME && unset PYTHONPATH && source /cvmfs/lhcb.cern.ch/lhcbdirac/lhcbdirac && python $PIDPERFSCRIPTSROOT/python/PIDPerfScripts/GetWGPFilesDirac.py %s %s %s %s'
        % (StripVer, MagPolarity, PartName, path)
    ]

    p = Popen(arguments, env=os.environ, stdout=PIPE, stderr=PIPE, shell=True)

    (cout, cerr) = p.communicate()
    print(cout)
    print(cerr)

    files = []

    with open(path + "/files.txt", 'rt') as fp:
        files = json.load(fp)

    #Remove the temporary file containing the files list
    os.remove(path + "/files.txt")
    os.rmdir

    return files


######################################################
### 2015 - 2018 ELECTRON ACCESS WITH THIS FUNCTION ###
######################################################
def GetElectronFiles(StripVer, MagPolarity, verbose=False):

    year = ''
    mag = ''

    if StripVer == 'Electron15':
        year = '2015'
    elif StripVer == 'Electron16':
        year = '2016'
    elif StripVer == 'Electron17':
        year = '2017'
    elif StripVer == 'Electron18':
        year = '2018'

    if MagPolarity == 'MagUp':
        mag = 'MU'
    elif MagPolarity == 'MagDown':
        mag = 'MD'

    newFileList = []

    if (year == "2015" or year == "2016"):
        suf = ""
        if (year == "2015"):
            suf = ".root"
        #Cuts applied to tag electron
        elif (year == "2016"):
            suf = "_TAGCUT.root"
        newFileList += [
            'root://eoslhcb.cern.ch//eos/lhcb/wg/PID/PIDCalib_' + year +
            '_electrons/pidcalib_BJpsiEE_' + mag + suf
        ]

    elif (year == "2017" or year == "2018"):
        newFileList += [
            'root://eoslhcb.cern.ch//eos/lhcb/wg/PID/PIDCalib_' + year +
            '_electrons/Turbo' + year + '_B2KJpsiEE_' + MagPolarity + '.root'
        ]

    return newFileList


#############################################
### RUN I FILE ACCESS WITH THIS FUNCTION ####
#############################################


def GetFiles(StripVer,
             MagPolarity,
             PartName,
             runMin=None,
             runMax=None,
             maxFiles=-1,
             verbose=False):
    #files = OverrideCalibDataStore.GetDictFiles(runMin,runMax,maxFiles,verbose)
    #if len(files) > 0:
    #   print "Using files from:" + os.getenv('OVERRIDECALIBDATASTORE')
    #   return files

    #else:
    #======================================================================
    # Create dictionary holding:
    # - Reconstruction version    ['RecoVer']
    # - np.array of:
    #        - MagUp run limits   ['UpRuns']
    #        - MagDown run limits ['DownRuns']
    #======================================================================
    files = []

    CheckStripVer(StripVer)

    CheckPartType(PartName)
    DataDict = GetRunDictionary(StripVer, PartName, verbose)

    #======================================================================
    # Determine min and max file indicies
    #======================================================================
    CheckMagPol(MagPolarity)
    IndexDict = GetMinMaxFileDictionary(DataDict, MagPolarity, runMin, runMax,
                                        maxFiles, verbose)

    from os import getenv

    RecoVer = GetRecoVer(StripVer)

    stv = StripVer
    if StripVer == '21_MCTuneV4':
        stv = '21'
    if StripVer == '21r1_MCTuneV4':
        stv = '21r1'
    if StripVer == '23_MCTuneV1':
        stv = '23'

    DataSetNameDict = GetDataSetNameDictionary(PartName)

    fname_protocol = ""
    fname_query = ""
    fname_extra = ""

    CalibDataProtocol = os.getenv("CALIBDATAURLPROTOCOL")
    CalibDataExtra = os.getenv("CALIBDATAEXTRA")

    # set the URL protocol (if applicable)
    if CalibDataProtocol is not None and CalibDataProtocol != "":
        fname_protocol = "{0}".format(CalibDataProtocol)

    if CalibDataExtra is not None and CalibDataExtra != "":
        fname_extra = "{0}".format(CalibDataExtra)

    vname_head = "CALIBDATASTORE"
    fname_head = os.getenv(vname_head)
    if fname_head is None:
        raise GetEnvError(
            "Cannot retrieve dataset, environmental variable %s has not been set."
            % vname_head)

    PartType = GetRealPartType(PartName)
    for i in range(IndexDict['minIndex'], IndexDict['maxIndex'] + 1):
        fname = ("{prtcol}//{extra}//{topdir}/Reco{reco}_DATA/{pol}/"
                 "{mother}_{part}_{pol}_Strip{strp}_{idx}.root").format(
                     prtcol=fname_protocol,
                     extra=fname_extra,
                     topdir=fname_head,
                     reco=RecoVer,
                     pol=MagPolarity,
                     mother=DataSetNameDict['MotherName'],
                     part=PartType,
                     strp=stv,
                     idx=i)
        files += [fname]
        print(files[-1])
    return files
