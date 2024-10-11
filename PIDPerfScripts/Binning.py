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
import ROOT
from PIDPerfScripts.Definitions import *

__all__ = ('CheckBinScheme', 'AddBinScheme', 'SetDefaultBinScheme',
           'AddUniformBins', 'AddBinBoundary', 'GetBinScheme')

_BinSchemes = {}
for trackType in GetPartTypes():
    for varName, varNameInDataSet in DataSetVariables().items():
        _BinSchemes.setdefault(trackType, {})
        _BinSchemes[trackType][varName] = {}


def CheckBinScheme(trackType, varName, schemeName, errorOnMissing=True):
    """Check whether the requested scheme name exists for the sepecified """
    """track type and bin variable.
If 'errorOnMissing' is True, then raises a KeyError exception upon finding """
    """a missing scheme name.
If 'errorOnMissing' is False, then returns False upon finding a missing """
    """scheme name.
Return True if the requested scheme name is found."""

    CheckPartType(trackType)
    CheckVarName(varName)
    schemeNames = list(_BinSchemes[trackType][varName].keys())
    if schemeName not in schemeNames:
        if errorOnMissing:
            schemeNames.sort()
            msg = ("Scheme name '{sname}' not in the list of bin schemes for "
                   "variable '{vname}'. Possible schemes are {snames}").format(
                       sname=schemeName,
                       vname=varName,
                       snames=str(schemeNames))
            raise KeyError(msg)
        else:
            return False
    return True


def AddBinScheme(trackType,
                 varName,
                 schemeName,
                 xMin,
                 xMax,
                 replaceCurrentScheme=False):
    """Adds a new scheme with the requested scheme name for the specified """
    """track type and bin variable with bin range [xMin, xMax].
If 'replaceCurrentScheme' is False, then a KeyError will be raised """
    """if a scheme with the same name already exists."""
    CheckPartType(trackType)
    CheckVarName(varName)
    schemeExists = CheckBinScheme(trackType, varName, schemeName, False)
    if schemeExists and not replaceCurrentScheme:
        msg = ("Scheme name '{sname}' already in the list of bin schemes "
               "for track type '{tname}', variable '{vname}'.").format(
                   sname=schemeName, tname=trackType, vname=varName)
        raise KeyError(msg)

    _BinSchemes[trackType][varName][schemeName] = ROOT.RooBinning(
        xMin, xMax, varName)


def SetDefaultBinScheme(trackType, varName, schemeName):
    """Set the default binning scheme for the specified track type and """
    """bin variable.
If 'muonUnBiased' is False, the default scheme for RICH calibration is set, """
    """otherwise the default scheme for muon calibration is set.
Raises a KeyError if a scheme with the requested name does not exist."""
    CheckBinScheme(trackType, varName, schemeName)
    _BinSchemes[trackType][varName]['default'] = ROOT.RooBinning(
        _BinSchemes[trackType][varName][schemeName], varName)


def AddUniformBins(trackType, varName, schemeName, nBins, xMin, xMax):
    """Adds 'nBins' bins, uniform in the range [xMin, xMax] to the binning """
    """scheme of the specified track type and bin variable.
Raises a KeyError if a scheme with the requested name does not exist."""
    CheckBinScheme(trackType, varName, schemeName)
    _BinSchemes[trackType][varName][schemeName].addUniform(nBins, xMin, xMax)


def AddBinBoundary(trackType, varName, schemeName, boundary):
    """Adds a new bin boundary to the binning scheme of the specified """
    """track type and bin variable.
Raises a KeyError if a scheme with the requested name does not exist."""
    CheckBinScheme(trackType, varName, schemeName)
    _BinSchemes[trackType][varName][schemeName].addBoundary(boundary)


def GetBinScheme(trackType, varName, schemeName=None):
    """Returns a copy of the requested binning scheme (we don't want the """
    """user to modify the original scheme).
If no scheme name is specified, then the default calibration scheme is """
    """used instead.
Raises a KeyError if a scheme with the requested name does not exist."""
    if schemeName is not None:
        CheckBinScheme(trackType, varName, schemeName)
        return ROOT.RooBinning(_BinSchemes[trackType][varName][schemeName],
                               varName)
    else:
        CheckBinScheme(trackType, varName, 'default')
        return ROOT.RooBinning(_BinSchemes[trackType][varName]['default'],
                               varName)


###########################################################################
######        Here, we make the default binning schemes              ######
###### The user can add more binning schemes using the above methods ######
###########################################################################

for trType in GetRICHPIDPartTypes() + GetMuonPIDPartTypes():
    # momentum
    AddBinScheme(trType, 'P', 'highres', 3000, 200000)
    AddUniformBins(trType, 'P', 'highres', 500, 3000, 200000)
    AddBinBoundary(trType, 'P', 'highres', 9300)  # R1 Kaon threshold
    AddBinBoundary(trType, 'P', 'highres', 15600)  # R2 Kaon threshold

    # momentum
    AddBinScheme(trType, 'Brunel_P', 'highres', 3000, 200000)
    AddUniformBins(trType, 'Brunel_P', 'highres', 500, 3000, 200000)
    AddBinBoundary(trType, 'Brunel_P', 'highres', 9300)  # R1 Kaon threshold
    AddBinBoundary(trType, 'Brunel_P', 'highres', 15600)  # R2 Kaon threshold

    # eta
    AddBinScheme(trType, 'ETA', 'highres', 1.5, 5)
    AddUniformBins(trType, 'ETA', 'highres', 500, 1.5, 5.0)

    # eta
    AddBinScheme(trType, 'Brunel_ETA', 'highres', 1.5, 5)
    AddUniformBins(trType, 'Brunel_ETA', 'highres', 500, 1.5, 5.0)

    # nTracks
    AddBinScheme(trType, 'nTracks', 'highres', 0, 500)
    AddUniformBins(trType, 'nTracks', 'highres', 500, 0, 500)

    # nTracks_Brunel (offline version)
    AddBinScheme(trType, 'nTracks_Brunel', 'highres', 0, 500)
    AddUniformBins(trType, 'nTracks_Brunel', 'highres', 500, 0, 500)

    AddBinScheme(trType, 'nSPDHits', 'highres', 0, 500)
    AddUniformBins(trType, 'nSPDHits', 'highres', 500, 0, 500)

    AddBinScheme(trType, 'nVeloClusters', 'highres', 0, 10000)
    AddUniformBins(trType, 'nVeloClusters', 'highres', 500, 0, 10000)

### DLL(K-pi), RICH (default schemes)

for trType in GetRICHPIDPartTypes():
    # momentum
    AddBinScheme(trType, 'P', 'DLLKpi', 3000, 100000)
    AddBinBoundary(trType, 'P', 'DLLKpi', 9300)  # R1 Kaon threshold
    AddBinBoundary(trType, 'P', 'DLLKpi', 15600)  # R2 Kaon threshold
    AddUniformBins(trType, 'P', 'DLLKpi', 15, 19000, 100000)

    # momentum
    AddBinScheme(trType, 'Brunel_P', 'DLLKpi', 3000, 100000)
    AddBinBoundary(trType, 'Brunel_P', 'DLLKpi', 9300)  # R1 Kaon threshold
    AddBinBoundary(trType, 'Brunel_P', 'DLLKpi', 15600)  # R2 Kaon threshold
    AddUniformBins(trType, 'Brunel_P', 'DLLKpi', 15, 19000, 100000)

    # eta
    AddBinScheme(trType, 'ETA', 'DLLKpi', 1.5, 5)
    AddUniformBins(trType, 'ETA', 'DLLKpi', 4, 1.5, 5)

    # eta
    AddBinScheme(trType, 'Brunel_ETA', 'DLLKpi', 1.5, 5)
    AddUniformBins(trType, 'Brunel_ETA', 'DLLKpi', 4, 1.5, 5)

    # nTracks
    AddBinScheme(trType, 'nTracks', 'DLLKpi', 0, 500)
    AddBinBoundary(trType, 'nTracks', 'DLLKpi', 50)
    AddBinBoundary(trType, 'nTracks', 'DLLKpi', 200)
    AddBinBoundary(trType, 'nTracks', 'DLLKpi', 300)

    # nTracks_Brunel (offline version)
    AddBinScheme(trType, 'nTracks_Brunel', 'DLLKpi', 0, 500)
    AddBinBoundary(trType, 'nTracks_Brunel', 'DLLKpi', 50)
    AddBinBoundary(trType, 'nTracks_Brunel', 'DLLKpi', 200)
    AddBinBoundary(trType, 'nTracks_Brunel', 'DLLKpi', 300)

    # nSPDHits

    AddBinScheme(trType, 'nSPDHits', 'DLLKpi', 0, 1000)
    AddUniformBins(trType, 'nSPDHits', 'DLLKpi', 5, 0, 1000)

    AddBinScheme(trType, 'nVeloClusters', 'DLLKpi', 0, 10000)
    AddUniformBins(trType, 'nVeloClusters', 'DLLKpi', 5, 0, 10000)

### DLL(K-pi), "MuonUnBiased" (default schemes)

for trType in GetMuonPIDPartTypes():
    # momentum
    AddBinScheme(trType, 'P', 'DLLKpi_MuonUnBiased', 3000, 100000)
    momBoundaries = (6000, 8000, 10000, 12000, 14500, 17500, 21500, 27000,
                     32000, 40000, 60000, 70000)
    for boundary in momBoundaries:
        AddBinBoundary(trType, 'P', 'DLLKpi_MuonUnBiased', boundary)

    # momentum
    AddBinScheme(trType, 'Brunel_P', 'DLLKpi_MuonUnBiased', 3000, 100000)
    momBoundaries = (6000, 8000, 10000, 12000, 14500, 17500, 21500, 27000,
                     32000, 40000, 60000, 70000)
    for boundary in momBoundaries:
        AddBinBoundary(trType, 'Brunel_P', 'DLLKpi_MuonUnBiased', boundary)

    # eta
    AddBinScheme(trType, 'ETA', 'DLLKpi_MuonUnBiased', 1.5, 5)
    AddUniformBins(trType, 'ETA', 'DLLKpi_MuonUnBiased', 4, 1.5, 5)

    # eta
    AddBinScheme(trType, 'Brunel_ETA', 'DLLKpi_MuonUnBiased', 1.5, 5)
    AddUniformBins(trType, 'Brunel_ETA', 'DLLKpi_MuonUnBiased', 4, 1.5, 5)

    # nTracks
    AddBinScheme(trType, 'nTracks', 'DLLKpi_MuonUnBiased', 0, 500)
    AddBinBoundary(trType, 'nTracks', 'DLLKpi_MuonUnBiased', 50)
    AddBinBoundary(trType, 'nTracks', 'DLLKpi_MuonUnBiased', 200)
    AddBinBoundary(trType, 'nTracks', 'DLLKpi_MuonUnBiased', 300)

    # nTracks_Brunel (offline version)
    AddBinScheme(trType, 'nTracks_Brunel', 'DLLKpi_MuonUnBiased', 0, 500)
    AddBinBoundary(trType, 'nTracks_Brunel', 'DLLKpi_MuonUnBiased', 50)
    AddBinBoundary(trType, 'nTracks_Brunel', 'DLLKpi_MuonUnBiased', 200)
    AddBinBoundary(trType, 'nTracks_Brunel', 'DLLKpi_MuonUnBiased', 300)

### RICH performance plots

for trType in ('K', 'Pi'):
    ## for K/pi ID/misID performance plots

    # momentum
    AddBinScheme(trType, 'P', 'PerfPlots_KPi', 2000, 100000)
    AddUniformBins(trType, 'P', 'PerfPlots_KPi', 40, 2000, 100000)

    # momentum
    AddBinScheme(trType, 'Brunel_P', 'PerfPlots_KPi', 2000, 100000)
    AddUniformBins(trType, 'Brunel_P', 'PerfPlots_KPi', 40, 2000, 100000)

    # eta
    AddBinScheme(trType, 'ETA', 'PerfPlots_KPi', 1.5, 5)
    AddUniformBins(trType, 'ETA', 'PerfPlots_KPi', 35, 1.5, 5)

    # eta
    AddBinScheme(trType, 'Brunel_ETA', 'PerfPlots_KPi', 1.5, 5)
    AddUniformBins(trType, 'Brunel_ETA', 'PerfPlots_KPi', 35, 1.5, 5)

    # nTracks
    AddBinScheme(trType, 'nTracks', 'PerfPlots_KPi', 0, 500)
    AddUniformBins(trType, 'nTracks', 'PerfPlots_KPi', 50, 0, 500)

    # nTracks_Brunel
    AddBinScheme(trType, 'nTracks_Brunel', 'PerfPlots_KPi', 0, 500)
    AddUniformBins(trType, 'nTracks_Brunel', 'PerfPlots_KPi', 50, 0, 500)

for trType in ('P', 'Pi'):
    ### for P/pi ID/misID performance plots

    # momentum
    AddBinScheme(trType, 'P', 'PerfPlots_PPi', 5000, 100000)
    AddUniformBins(trType, 'P', 'PerfPlots_PPi', 38, 5000, 100000)

    # momentum
    AddBinScheme(trType, 'Brunel_P', 'PerfPlots_PPi', 5000, 100000)
    AddUniformBins(trType, 'Brunel_P', 'PerfPlots_PPi', 38, 5000, 100000)

    # eta
    AddBinScheme(trType, 'ETA', 'PerfPlots_PPi', 1.5, 5)
    AddUniformBins(trType, 'ETA', 'PerfPlots_PPi', 35, 1.5, 5)

    # eta
    AddBinScheme(trType, 'Brunel_ETA', 'PerfPlots_PPi', 1.5, 5)
    AddUniformBins(trType, 'Brunel_ETA', 'PerfPlots_PPi', 35, 1.5, 5)

    # nTracks
    AddBinScheme(trType, 'nTracks', 'PerfPlots_PPi', 0, 500)
    AddUniformBins(trType, 'nTracks', 'PerfPlots_PPi', 50, 0, 500)

    # nTracks_Brunel
    AddBinScheme(trType, 'nTracks_Brunel', 'PerfPlots_PPi', 0, 500)
    AddUniformBins(trType, 'nTracks_Brunel', 'PerfPlots_PPi', 50, 0, 500)

for trType in ('Mu', 'K_MuonUnBiased'):
    ### for Mu/K ID/misID performance plots

    # momentum
    AddBinScheme(trType, 'P', 'PerfPlots_MuK_MuonUnBiased', 2000, 100000)
    AddUniformBins(trType, 'P', 'PerfPlots_MuK_MuonUnBiased', 40, 2000, 100000)

    # momentum
    AddBinScheme(trType, 'Brunel_P', 'PerfPlots_MuK_MuonUnBiased', 2000,
                 100000)
    AddUniformBins(trType, 'Brunel_P', 'PerfPlots_MuK_MuonUnBiased', 40, 2000,
                   100000)

    # eta
    AddBinScheme(trType, 'ETA', 'PerfPlots_MuK_MuonUnBiased', 1.5, 5)
    AddUniformBins(trType, 'ETA', 'PerfPlots_MuK_MuonUnBiased', 35, 1.5, 5)

    # eta
    AddBinScheme(trType, 'Brunel_ETA', 'PerfPlots_MuK_MuonUnBiased', 1.5, 5)
    AddUniformBins(trType, 'Brunel_ETA', 'PerfPlots_MuK_MuonUnBiased', 35, 1.5,
                   5)

    # nTracks
    AddBinScheme(trType, 'nTracks', 'PerfPlots_MuK_MuonUnBiased', 0, 500)
    AddUniformBins(trType, 'nTracks', 'PerfPlots_MuK_MuonUnBiased', 50, 0, 500)

    # nTracks_Brunel
    AddBinScheme(trType, 'nTracks_Brunel', 'PerfPlots_MuK_MuonUnBiased', 0,
                 500)
    AddUniformBins(trType, 'nTracks_Brunel', 'PerfPlots_MuK_MuonUnBiased', 50,
                   0, 500)

for trType in ('Mu', 'Pi_MuonUnBiased'):
    ### for Mu/pi ID/misID performance plots

    # momentum
    AddBinScheme(trType, 'P', 'PerfPlots_MuPi_MuonUnBiased', 2000, 100000)
    AddUniformBins(trType, 'P', 'PerfPlots_MuPi_MuonUnBiased', 40, 2000,
                   100000)

    # momentum
    AddBinScheme(trType, 'Brunel_P', 'PerfPlots_MuPi_MuonUnBiased', 2000,
                 100000)
    AddUniformBins(trType, 'Brunel_P', 'PerfPlots_MuPi_MuonUnBiased', 40, 2000,
                   100000)

    # eta
    AddBinScheme(trType, 'ETA', 'PerfPlots_MuPi_MuonUnBiased', 1.5, 5)
    AddUniformBins(trType, 'ETA', 'PerfPlots_MuPi_MuonUnBiased', 35, 1.5, 5)

    # eta
    AddBinScheme(trType, 'Brunel_ETA', 'PerfPlots_MuPi_MuonUnBiased', 1.5, 5)
    AddUniformBins(trType, 'Brunel_ETA', 'PerfPlots_MuPi_MuonUnBiased', 35,
                   1.5, 5)

    # nTracks
    AddBinScheme(trType, 'nTracks', 'PerfPlots_MuPi_MuonUnBiased', 0, 500)
    AddUniformBins(trType, 'nTracks', 'PerfPlots_MuPi_MuonUnBiased', 50, 0,
                   500)

    # nTracks_Brunel
    AddBinScheme(trType, 'nTracks_Brunel', 'PerfPlots_MuPi_MuonUnBiased', 0,
                 500)
    AddUniformBins(trType, 'nTracks_Brunel', 'PerfPlots_MuPi_MuonUnBiased', 50,
                   0, 500)

for trType in ('Mu', 'P_MuonUnBiased'):
    ### for Mu/pi ID/misID performance plots

    # momentum
    AddBinScheme(trType, 'P', 'PerfPlots_MuP_MuonUnBiased', 2000, 100000)
    AddUniformBins(trType, 'P', 'PerfPlots_MuP_MuonUnBiased', 40, 2000, 100000)

    # momentum
    AddBinScheme(trType, 'Brunel_P', 'PerfPlots_MuP_MuonUnBiased', 2000,
                 100000)
    AddUniformBins(trType, 'Brunel_P', 'PerfPlots_MuP_MuonUnBiased', 40, 2000,
                   100000)

    # eta
    AddBinScheme(trType, 'ETA', 'PerfPlots_MuP_MuonUnBiased', 1.5, 5)
    AddUniformBins(trType, 'ETA', 'PerfPlots_MuP_MuonUnBiased', 35, 1.5, 5)

    # eta
    AddBinScheme(trType, 'Brunel_ETA', 'PerfPlots_MuP_MuonUnBiased', 1.5, 5)
    AddUniformBins(trType, 'Brunel_ETA', 'PerfPlots_MuP_MuonUnBiased', 35, 1.5,
                   5)

    # nTracks
    AddBinScheme(trType, 'nTracks', 'PerfPlots_MuP_MuonUnBiased', 0, 500)
    AddUniformBins(trType, 'nTracks', 'PerfPlots_MuP_MuonUnBiased', 50, 0, 500)

    # nTracks_Brunel
    AddBinScheme(trType, 'nTracks_Brunel', 'PerfPlots_MuP_MuonUnBiased', 0,
                 500)
    AddUniformBins(trType, 'nTracks_Brunel', 'PerfPlots_MuP_MuonUnBiased', 50,
                   0, 500)

for trType in ('e', 'Pi'):
    ### for e/pi ID/misID performance plots

    # momentum
    AddBinScheme(trType, 'P', 'PerfPlots_ePi', 5000, 100000)
    AddUniformBins(trType, 'P', 'PerfPlots_ePi', 38, 5000, 100000)

    # momentum
    AddBinScheme(trType, 'Brunel_P', 'PerfPlots_ePi', 5000, 100000)
    AddUniformBins(trType, 'Brunel_P', 'PerfPlots_ePi', 38, 5000, 100000)

    # eta
    AddBinScheme(trType, 'ETA', 'PerfPlots_ePi', 1.5, 5)
    AddUniformBins(trType, 'ETA', 'PerfPlots_ePi', 35, 1.5, 5)

    # eta
    AddBinScheme(trType, 'Brunel_ETA', 'PerfPlots_ePi', 1.5, 5)
    AddUniformBins(trType, 'Brunel_ETA', 'PerfPlots_ePi', 35, 1.5, 5)

    # nTracks
    AddBinScheme(trType, 'nTracks', 'PerfPlots_ePi', 0, 500)
    AddUniformBins(trType, 'nTracks', 'PerfPlots_ePi', 50, 0, 500)

    # nTracks_Brunel
    AddBinScheme(trType, 'nTracks_Brunel', 'PerfPlots_ePi', 0, 500)
    AddUniformBins(trType, 'nTracks_Brunel', 'PerfPlots_ePi', 50, 0, 500)

for trType in GetPartTypes():
    AddBinScheme(trType, 'P', 'PerfPlots_%s' % (trType), 0, 100000)
    AddUniformBins(trType, 'P', 'PerfPlots_%s' % (trType), 20, 0, 100000)

    AddBinScheme(trType, 'Brunel_P', 'PerfPlots_%s' % (trType), 0, 100000)
    AddUniformBins(trType, 'Brunel_P', 'PerfPlots_%s' % (trType), 20, 0,
                   100000)

    AddBinScheme(trType, 'PT', 'PerfPlots_%s' % (trType), 0, 15000)
    AddUniformBins(trType, 'PT', 'PerfPlots_%s' % (trType), 20, 0, 15000)

    AddBinScheme(trType, 'Brunel_PT', 'PerfPlots_%s' % (trType), 0, 15000)
    AddUniformBins(trType, 'Brunel_PT', 'PerfPlots_%s' % (trType), 20, 0,
                   15000)

    AddBinScheme(trType, 'ETA', 'PerfPlots_%s' % (trType), 1.5, 5.0)
    AddUniformBins(trType, 'ETA', 'PerfPlots_%s' % (trType), 20, 1.5, 5.0)

    AddBinScheme(trType, 'Brunel_ETA', 'PerfPlots_%s' % (trType), 1.5, 5.0)
    AddUniformBins(trType, 'Brunel_ETA', 'PerfPlots_%s' % (trType), 20, 1.5,
                   5.0)

    AddBinScheme(trType, 'PHI', 'PerfPlots_%s' % (trType), -3.14159, 3.14159)
    AddUniformBins(trType, 'PHI', 'PerfPlots_%s' % (trType), 20, -3.14159,
                   3.14159)

    AddBinScheme(trType, 'runNumber', 'PerfPlots_%s' % (trType), 87660, 104300)
    AddUniformBins(trType, 'runNumber', 'PerfPlots_%s' % (trType), 100, 87660,
                   104300)

    for realPart in [p.lower() for p in GetRealPartTypes()]:
        if realPart == "k":
            realPart = "K"

        if realPart != "pi":
            AddBinScheme(trType, 'DLL%s' % (realPart),
                         'PerfPlots_%s' % (trType), -100, 100)
            AddUniformBins(trType, 'DLL%s' % (realPart),
                           'PerfPlots_%s' % (trType), 100, -100, 100)

        AddBinScheme(trType, 'MC12TuneV2_ProbNN%s' % (realPart),
                     'PerfPlots_%s' % (trType), 0, 1)
        AddUniformBins(trType, 'MC12TuneV2_ProbNN%s' % (realPart),
                       'PerfPlots_%s' % (trType), 20, 0, 1)

        AddBinScheme(trType, 'Brunel_MC12TuneV2_ProbNN%s' % (realPart),
                     'PerfPlots_%s' % (trType), 0, 1)
        AddUniformBins(trType, 'Brunel_MC12TuneV2_ProbNN%s' % (realPart),
                       'PerfPlots_%s' % (trType), 20, 0, 1)

        AddBinScheme(trType, 'MC12TuneV3_ProbNN%s' % (realPart),
                     'PerfPlots_%s' % (trType), 0, 1)
        AddUniformBins(trType, 'MC12TuneV3_ProbNN%s' % (realPart),
                       'PerfPlots_%s' % (trType), 20, 0, 1)

        AddBinScheme(trType, 'Brunel_MC12TuneV3_ProbNN%s' % (realPart),
                     'PerfPlots_%s' % (trType), 0, 1)
        AddUniformBins(trType, 'Brunel_MC12TuneV3_ProbNN%s' % (realPart),
                       'PerfPlots_%s' % (trType), 20, 0, 1)

        AddBinScheme(trType, 'MC12TuneV4_ProbNN%s' % (realPart),
                     'PerfPlots_%s' % (trType), 0, 1)
        AddUniformBins(trType, 'MC12TuneV4_ProbNN%s' % (realPart),
                       'PerfPlots_%s' % (trType), 20, 0, 1)

        AddBinScheme(trType, 'Brunel_MC12TuneV4_ProbNN%s' % (realPart),
                     'PerfPlots_%s' % (trType), 0, 1)
        AddUniformBins(trType, 'Brunel_MC12TuneV4_ProbNN%s' % (realPart),
                       'PerfPlots_%s' % (trType), 20, 0, 1)

        AddBinScheme(trType, 'MC15TuneV1_ProbNN%s' % (realPart),
                     'PerfPlots_%s' % (trType), 0, 1)
        AddUniformBins(trType, 'MC15TuneV1_ProbNN%s' % (realPart),
                       'PerfPlots_%s' % (trType), 20, 0, 1)

        AddBinScheme(trType, 'Brunel_MC15TuneV1_ProbNN%s' % (realPart),
                     'PerfPlots_%s' % (trType), 0, 1)
        AddUniformBins(trType, 'Brunel_MC15TuneV1_ProbNN%s' % (realPart),
                       'PerfPlots_%s' % (trType), 20, 0, 1)

        AddBinScheme(trType, 'MC15TuneDNNV1_ProbNN%s' % (realPart),
                     'PerfPlots_%s' % (trType), 0, 1)
        AddUniformBins(trType, 'MC15TuneDNNV1_ProbNN%s' % (realPart),
                       'PerfPlots_%s' % (trType), 20, 0, 1)

        AddBinScheme(trType, 'Brunel_MC15TuneDNNV1_ProbNN%s' % (realPart),
                     'PerfPlots_%s' % (trType), 0, 1)
        AddUniformBins(trType, 'Brunel_MC15TuneDNNV1_ProbNN%s' % (realPart),
                       'PerfPlots_%s' % (trType), 20, 0, 1)

        AddBinScheme(trType, 'MC15TuneFLAT4dV1_ProbNN%s' % (realPart),
                     'PerfPlots_%s' % (trType), 0, 1)
        AddUniformBins(trType, 'MC15TuneFLAT4dV1_ProbNN%s' % (realPart),
                       'PerfPlots_%s' % (trType), 20, 0, 1)

        AddBinScheme(trType, 'Brunel_MC15TuneFLAT4dV1_ProbNN%s' % (realPart),
                     'PerfPlots_%s' % (trType), 0, 1)
        AddUniformBins(trType, 'Brunel_MC15TuneFLAT4dV1_ProbNN%s' % (realPart),
                       'PerfPlots_%s' % (trType), 20, 0, 1)

        AddBinScheme(trType, 'MC15TuneCatBoostV1_ProbNN%s' % (realPart),
                     'PerfPlots_%s' % (trType), 0, 1)
        AddUniformBins(trType, 'MC15TuneCatBoostV1_ProbNN%s' % (realPart),
                       'PerfPlots_%s' % (trType), 20, 0, 1)

        AddBinScheme(trType, 'Brunel_MC15TuneCatBoostV1_ProbNN%s' % (realPart),
                     'PerfPlots_%s' % (trType), 0, 1)
        AddUniformBins(trType,
                       'Brunel_MC15TuneCatBoostV1_ProbNN%s' % (realPart),
                       'PerfPlots_%s' % (trType), 20, 0, 1)

        if realPart != "mu":
            AddBinScheme(trType, 'RICHThreshold_%s' % (realPart),
                         'PerfPlots_%s' % (trType), -2, 2)
            AddUniformBins(trType, 'RICHThreshold_%s' % (realPart),
                           'PerfPlots_%s' % (trType), 4, -2, 2)

    #Deuteron PID (offline only)
    AddBinScheme(trType, 'Brunel_MC15TuneV1_ProbNNd',
                 'PerfPlots_%s' % (trType), 0, 1)
    AddUniformBins(trType, 'Brunel_MC15TuneV1_ProbNNd',
                   'PerfPlots_%s' % (trType), 20, 0, 1)

    AddBinScheme(
        trType,
        'DLLe',
        'PerfPlots_%s' % (trType),
        -15,
        20,
        replaceCurrentScheme=True)
    AddUniformBins(trType, 'DLLe', 'PerfPlots_%s' % (trType), 100, -15, 20)

    AddBinScheme(
        trType,
        'Brunel_DLLe',
        'PerfPlots_%s' % (trType),
        -15,
        20,
        replaceCurrentScheme=True)
    AddUniformBins(trType, 'Brunel_DLLe', 'PerfPlots_%s' % (trType), 100, -15,
                   20)

    AddBinScheme(trType, 'InMuonAcc', 'PerfPlots_%s' % (trType), 0, 2)
    AddUniformBins(trType, 'InMuonAcc', 'PerfPlots_%s' % (trType), 2, 0, 2)

    AddBinScheme(trType, 'Brunel_InMuonAcc', 'PerfPlots_%s' % (trType), 0, 2)
    AddUniformBins(trType, 'Brunel_InMuonAcc', 'PerfPlots_%s' % (trType), 2, 0,
                   2)

    AddBinScheme(trType, 'IsMuon', 'PerfPlots_%s' % (trType), 0, 2)
    AddUniformBins(trType, 'IsMuon', 'PerfPlots_%s' % (trType), 2, 0, 2)

    AddBinScheme(trType, 'Brunel_IsMuon', 'PerfPlots_%s' % (trType), 0, 2)
    AddUniformBins(trType, 'Brunel_IsMuon', 'PerfPlots_%s' % (trType), 2, 0, 2)

    AddBinScheme(trType, 'IsMuonLoose', 'PerfPlots_%s' % (trType), 0, 2)
    AddUniformBins(trType, 'IsMuonLoose', 'PerfPlots_%s' % (trType), 2, 0, 2)

    AddBinScheme(trType, 'Brunel_IsMuonLoose', 'PerfPlots_%s' % (trType), 0, 2)
    AddUniformBins(trType, 'Brunel_IsMuonLoose', 'PerfPlots_%s' % (trType), 2,
                   0, 2)

    AddBinScheme(trType, 'nShared', 'PerfPlots_%s' % (trType), 0, 200)
    AddUniformBins(trType, 'nShared', 'PerfPlots_%s' % (trType), 50, 0, 200)

    AddBinScheme(trType, 'Brunel_nShared', 'PerfPlots_%s' % (trType), 0, 200)
    AddUniformBins(trType, 'Brunel_nShared', 'PerfPlots_%s' % (trType), 50, 0,
                   200)

    AddBinScheme(trType, 'RICHAerogelUsed', 'PerfPlots_%s' % (trType), 0, 2)
    AddUniformBins(trType, 'RICHAerogelUsed', 'PerfPlots_%s' % (trType), 2, 0,
                   2)

    AddBinScheme(trType, 'Brunel_RICHAerogelUsed', 'PerfPlots_%s' % (trType),
                 0, 2)
    AddUniformBins(trType, 'Brunel_RICHAerogelUsed', 'PerfPlots_%s' % (trType),
                   2, 0, 2)

    AddBinScheme(trType, 'RICH1GasUsed', 'PerfPlots_%s' % (trType), 0, 2)
    AddUniformBins(trType, 'RICH1GasUsed', 'PerfPlots_%s' % (trType), 2, 0, 2)

    AddBinScheme(trType, 'Brunel_RICH1GasUsed', 'PerfPlots_%s' % (trType), 0,
                 2)
    AddUniformBins(trType, 'Brunel_RICH1GasUsed', 'PerfPlots_%s' % (trType), 2,
                   0, 2)

    AddBinScheme(trType, 'RICH2GasUsed', 'PerfPlots_%s' % (trType), 0, 2)
    AddUniformBins(trType, 'RICH2GasUsed', 'PerfPlots_%s' % (trType), 2, 0, 2)

    AddBinScheme(trType, 'Brunel_RICH2GasUsed', 'PerfPlots_%s' % (trType), 0,
                 2)
    AddUniformBins(trType, 'Brunel_RICH2GasUsed', 'PerfPlots_%s' % (trType), 2,
                   0, 2)

    AddBinScheme(trType, 'HasRich', 'PerfPlots_%s' % (trType), 0, 2)
    AddUniformBins(trType, 'HasRich', 'PerfPlots_%s' % (trType), 2, 0, 2)

    AddBinScheme(trType, 'Brunel_HasRich', 'PerfPlots_%s' % (trType), 0, 2)
    AddUniformBins(trType, 'Brunel_HasRich', 'PerfPlots_%s' % (trType), 2, 0,
                   2)

    AddBinScheme(trType, 'HasCalo', 'PerfPlots_%s' % (trType), 0, 2)
    AddUniformBins(trType, 'HasCalo', 'PerfPlots_%s' % (trType), 2, 0, 2)

    AddBinScheme(trType, 'Brunel_HasCalo', 'PerfPlots_%s' % (trType), 0, 2)
    AddUniformBins(trType, 'Brunel_HasCalo', 'PerfPlots_%s' % (trType), 2, 0,
                   2)

    #AddBinScheme(trType, 'HasBremAdded', 'PerfPlots_%s'%(trType), 0, 2)
    #AddUniformBins(trType, 'HasBremAdded', 'PerfPlots_%s'%(trType), 2, 0, 2)

    #AddBinScheme(trType, 'CaloRegion', 'PerfPlots_%s'%(trType), 2, 5)
    #AddUniformBins(trType, 'CaloRegion', 'PerfPlots_%s'%(trType), 3, 2, 5)

    AddBinScheme(trType, 'nTracks', 'PerfPlots_%s' % (trType), 0, 500)
    AddUniformBins(trType, 'nTracks', 'PerfPlots_%s' % (trType), 20, 0, 500)

    AddBinScheme(trType, 'nTracks_Brunel', 'PerfPlots_%s' % (trType), 0, 500)
    AddUniformBins(trType, 'nTracks_Brunel', 'PerfPlots_%s' % (trType), 20, 0,
                   500)

    AddBinScheme(trType, 'nSPDHits', 'PerfPlots_%s' % (trType), 0, 1000)
    AddUniformBins(trType, 'nSPDHits', 'PerfPlots_%s' % (trType), 20, 0, 1000)

    AddBinScheme(trType, 'nVeloClusters', 'PerfPlots_%s' % (trType), 0, 10000)
    AddUniformBins(trType, 'nVeloClusters', 'PerfPlots_%s' % (trType), 20, 0,
                   10000)

    AddBinScheme(trType, 'nRich1Hits', 'PerfPlots_%s' % (trType), 0, 10000)
    AddUniformBins(trType, 'nRich1Hits', 'PerfPlots_%s' % (trType), 20, 0,
                   10000)

    AddBinScheme(trType, 'nRich2Hits', 'PerfPlots_%s' % (trType), 0, 8000)
    AddUniformBins(trType, 'nRich2Hits', 'PerfPlots_%s' % (trType), 20, 0,
                   8000)

### set the default binning schemes
for trType in GetRICHPIDPartTypes():
    SetDefaultBinScheme(trType, 'P', 'DLLKpi')
    SetDefaultBinScheme(trType, 'Brunel_P', 'DLLKpi')
    SetDefaultBinScheme(trType, 'ETA', 'DLLKpi')
    SetDefaultBinScheme(trType, 'Brunel_ETA', 'DLLKpi')
    SetDefaultBinScheme(trType, 'nTracks', 'DLLKpi')
    SetDefaultBinScheme(trType, 'nTracks_Brunel', 'DLLKpi')
    SetDefaultBinScheme(trType, 'nSPDHits', 'DLLKpi')
    SetDefaultBinScheme(trType, 'nVeloClusters', 'DLLKpi')

for trType in GetMuonPIDPartTypes():
    SetDefaultBinScheme(trType, 'P', 'DLLKpi_MuonUnBiased')
    SetDefaultBinScheme(trType, 'Brunel_P', 'DLLKpi_MuonUnBiased')
    SetDefaultBinScheme(trType, 'ETA', 'DLLKpi_MuonUnBiased')
    SetDefaultBinScheme(trType, 'Brunel_ETA', 'DLLKpi_MuonUnBiased')
    SetDefaultBinScheme(trType, 'nTracks', 'DLLKpi_MuonUnBiased')
    SetDefaultBinScheme(trType, 'nTracks_Brunel', 'DLLKpi_MuonUnBiased')

#for trType in GetProtonPIDPartTypes():
#    SetDefaultBinScheme(trType, 'P', 'DLLKpi')
#    SetDefaultBinScheme(trType, 'ETA', 'DLLKpi')
#    SetDefaultBinScheme(trType, 'nTracks', 'DLLKpi')
