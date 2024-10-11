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
from __future__ import absolute_import
import re
import math

from PIDPerfScripts.Exceptions import *
from .TupleDataset import *


def BannedVariables(var=None):
    h = {
        'V2ProbNNK': '{particle}_V2ProbNNK',
        'V2ProbNNpi': '{particle}_V2ProbNNpi',
        'V2ProbNNp': '{particle}_V2ProbNNp',
        'V2ProbNNmu': '{particle}_V2ProbNNmu',
        'V2ProbNNe': '{particle}_V2ProbNNe',
        'V2ProbNNghost': '{particle}_V2ProbNNghost',
        'V3ProbNNK': '{particle}_V3ProbNNK',
        'V3ProbNNpi': '{particle}_V3ProbNNpi',
        'V3ProbNNp': '{particle}_V3ProbNNp',
        'V3ProbNNmu': '{particle}_V3ProbNNmu',
        'V3ProbNNe': '{particle}_V3ProbNNe',
        'V3ProbNNghost': '{particle}_V3ProbNNghost'
    }
    if var == None:
        return h
    else:
        return h[var]


#Load the DataSetVariables from the list defined in TupleDataset
def DataSetVariables(var=None):
    h = vars_dataset
    if var == None:
        return h
    else:
        return h[var]


def GetVarNames():
    varArr = [
        varName for varName, varNameInDataset in DataSetVariables().items()
    ]
    varArr.sort()
    return varArr


def GetRICHPIDRealPartTypes():
    return ("K", "Pi", "P", "e")


def GetRICHPIDPartTypes():
    return GetRICHPIDRealPartTypes() + GetProtonPIDPartTypes(
    ) + GetDDPIDPartTypes() + ("K_DsPhi", "K_DSt3Pi", "K_Phi", "K_Omega",
                               "Pi_KS", "Pi_DSt3Pi", "e_B_Jpsi")


def GetMuonPIDRealPartTypes():
    return ("Mu", )


def GetProtonPIDPartTypes():
    return ("P_LcfB", "P_IncLc", "P_TotLc") + ("P_Sigmac0", "P_Sigmacpp",
                                               "P_B_Jpsi")


def GetProtonPIDRealPartTypes():
    return ("P", )


def GetDDPIDPartTypes():
    return ("Pi_DD", "K_DD", "P_DD")


def GetMuonPIDPartTypes():
    return GetMuonPIDRealPartTypes() + tuple(
        ['{0}_MuonUnBiased'.format(p)
         for p in GetRICHPIDRealPartTypes()]) + ("Mu_B_Jpsi", "Mu_nopt")


def GetRealPartTypes():
    return GetRICHPIDRealPartTypes() + GetMuonPIDRealPartTypes()


def GetPartTypes():
    return GetRICHPIDPartTypes() + GetMuonPIDPartTypes()


def GetRealPartType(PartName):
    CheckPartType(PartName)
    if PartName == 'K' or PartName == 'K_MuonUnBiased' or PartName == 'K_DsPhi' or PartName == 'K_DSt3Pi' or PartName == 'K_Phi' or PartName == 'K_Omega' or PartName == 'K_DD':
        return 'K'
    elif PartName == 'Pi' or PartName == 'Pi_MuonUnBiased' or PartName == 'Pi_KS' or PartName == 'Pi_DSt3Pi' or PartName == 'Pi_DD':
        return 'Pi'
    elif PartName == 'P' or PartName == 'P_MuonUnBiased' or PartName == 'P_LcfB' or PartName == 'P_IncLc' or PartName == 'P_TotLc' or PartName == 'P_Sigmac0' or PartName == 'P_Sigmacpp' or PartName == 'P_B_Jpsi' or PartName == 'P_DD':
        return 'P'
    elif PartName == 'e' or PartName == 'e_MuonUnBiased' or PartName == 'e_B_Jpsi':
        return 'e'
    else:
        return 'Mu'


# N.B. The following method will need to be changed once the
# Lambda_c proton samples are included
def GetMotherName(PartName):
    #PartType = GetRealPartType(PartName)
    if PartName == 'K' or PartName == 'Pi':
        return 'DSt'
    elif PartName == 'K_DSt3Pi' or PartName == 'Pi_DSt3Pi':
        return 'DSt3Pi'
    elif PartName == 'K_DsPhi':
        return 'DsPhi'
    elif PartName == 'K_Phi':
        return 'Phi'
    elif PartName == 'K_Omega':
        return 'Omega'
    elif PartName == 'K_DD':
        return 'Omega_DDD'
    elif PartName == 'K_KS' or PartName == 'Pi_KS':
        return 'KS'
    elif PartName == 'Pi_DD':
        return 'KS_DD'
    elif PartName == 'P':
        return 'Lam0'
    elif PartName == 'P_DD':
        return 'Lam0_DD'
    elif PartName == 'P_Sigmac0':
        return 'Sigmac0'
    elif PartName == 'P_Sigmacpp':
        return 'Sigmacpp'
    elif PartName == 'P_B_Jpsi' or PartName == 'e_B_Jpsi' or PartName == 'Mu_B_Jpsi':
        return 'B_Jpsi'
    elif PartName == 'e':
        return 'Jpsi'
    elif PartName == 'K_MuonUnBiased' or PartName == 'Pi_MuonUnBiased':
        return 'DSt_MuonUnBiased'
    elif PartName == 'P_MuonUnBiased':
        return 'Lam0_MuonUnBiased'
    elif PartName == 'e_MuonUnBiased':
        return 'Jpsi'
    elif PartName == 'P_LcfB':
        return 'LcfB'
    elif PartName == 'P_TotLc':
        return 'TotLc'
    elif PartName == 'P_IncLc':
        return 'IncLc'
    #Jpsi muon sample with no pT cut
    elif PartName == 'Mu_nopt':
        return 'Jpsi_nopt'
    else:
        return 'Jpsi'


# N.B. The following method will need to be changed once the
# Lambda_c proton samples are included - its going to come from the mother but not for electrons and muons as the are both from jpsis
def GetWorkspaceName(PartName):
    PartType = GetRealPartType(PartName)
    MType = GetMotherName(PartName)
    if MType == 'DSt' or MType == 'Dst_MuonUnBiased':
        return 'RSDStCalib'
    elif MType == 'Lam0' or MType == 'Lam0_MuonUnBiased':
        return 'Lam0Calib'
    elif PartType == 'e':
        return 'JpsieeCalib'
    elif MType == 'LcfB':
        return 'IncLcCalib'
    elif MType == 'TotLc':
        return 'IncLcCalib'
    elif MType == 'IncLc':
        return 'IncLcCalib'
    elif MType == 'Jpsi':
        return 'JpsiCalib'
    else:
        print('No workspace found,this is working group production')


def IsMuonUnBiased(PartName):
    CheckPartType(PartName)
    if PartName in GetMuonPIDPartTypes():
        return True
    else:
        return False


def GetRecoVer(StripVer):
    CheckStripVer(StripVer)
    if StripVer == '13b':
        return 10
    elif StripVer == '15':
        return 11
    elif StripVer == '17':
        return 12
    elif StripVer == '22':
        return '15a'
    elif StripVer == '23Val':
        return '15a'
    elif StripVer == '23':
        return '15a'
    elif StripVer == '5TeV':
        return '15a'
    elif StripVer == '21_MCTuneV4':
        return '14V4'
    elif StripVer == '21r1_MCTuneV4':
        return '14V4'
    elif StripVer == '23_MCTuneV1':
        return '15aV1'
    elif StripVer == '26':
        return '16'
    else:
        return 14


def GetKnownFunctions():
    ret = dir(math)
    ret += ["abs", "fabs"]
    return ret


def CheckCutVarsInTree(cuts, tree):
    if cuts in [None, ""]:
        return True
    math_chars = re.compile("\\s*([-=|&\\(\\)\\+\\*/<>\\[\\]!\\s]+)\\s*")
    split_cuts = re.split(math_chars, cuts)
    #print "CheckCutVarsInTree():", split_cuts
    for var in split_cuts:
        found = False
        if re.match(math_chars, var):
            #print var,"is mathematical"
            continue
        elif var in GetKnownFunctions():
            #print var,"is a function"
            continue
        else:
            try:
                float(var)
                #print var,"is a number"
                continue
            except ValueError:
                pass
            #print var,"is not mathematical or a function"
            tree_vars = tree.GetListOfLeaves()
            for i in range(len(tree_vars)):
                if tree_vars[i].GetName() == var:
                    #print var,"==",tree_vars[i].GetName()
                    found = True
                    break
        if not found:
            raise ValueError(
                "Variable '%s' in cut string '%s' not found in TTree %s!" %
                (var, cuts, tree.GetName()))
            return False
    #print "Everything is bon!"
    return True


def CheckCuts(cuts, strip):
    simple_cuts = cuts
    for i in "()[]&!|><=:+-*/":
        simple_cuts = simple_cuts.replace(i, " ")

    simple_cuts = [x for x in simple_cuts.split(" ") if x != ""]
    valid_varibles = [x for x, y in DataSetVariables().items()]
    banned_varibles = [x for x, y in BannedVariables().items()]

    valid_varibles.sort()

    for var in simple_cuts:
        try:
            float(var)
        except ValueError:
            if var not in valid_varibles:
                print("'%s' is not a valid variable" % var)
                print("Known variables are:")
                print(valid_varibles)
                return False
            if strip == '23_MCTuneV1' or strip == '26':
                if var in banned_varibles:
                    print("")
                    print(
                        "**********************************************************"
                    )
                    print("'%s' is banned for RUN2" % var)
                    print(
                        "Please use either MC15TuneV1 (V1ProbNN) or MC12TuneV4 (V4ProbNN)"
                    )
                    print("For futher information, please check ")
                    print(
                        "https://twiki.cern.ch/twiki/bin/view/LHCb/PIDCalibPackage#PIDCalib_sample_for_2015_and_201"
                    )
                    print(
                        "*********************************************************"
                    )
                    print("")
                    print("")
                    return False
    return True


def FlatternPlots(Plots):
    ret = []
    if isinstance(Plots, list):
        for plot in Plots:
            ret.extend(FlatternPlots(plot))
    else:
        return [Plots]
    return ret


def CheckRealPartType(PartName):
    ValidPartNames = GetRealPartTypes()
    if PartName not in ValidPartNames:
        msg = ("Invalid particle type '{0}'. "
               "Allowed types are {1}").format(PartName, str(ValidPartNames))
        raise TypeError(msg)


def CheckPartType(PartName):
    ValidPartNames = GetPartTypes()
    if PartName not in ValidPartNames:
        msg = ("Invalid particle type '{0}'. "
               "Allowed types are {1}").format(PartName, str(ValidPartNames))
        raise TypeError(msg)


def CheckMagPol(MagPol):
    ValidMagPols = ("MagUp", "MagDown")
    if MagPol not in ValidMagPols:
        msg = ("Invalid magnet polarity '{0}'. "
               "Allowed polarities are {1}").format(MagPol, str(ValidMagPols))
        raise TypeError(msg)


def CheckStripVer(StripVer):
    ValidStripVers = ("13b", "15", "17", "20", "20r1", "20_MCTuneV2",
                      "23_MCTuneV1", "26", "20r1_MCTuneV2", "20_MCTunev3",
                      "21r1", "21", "22", "23Val", "23", "5TeV", "21_MCTuneV4",
                      "21r1_MCTuneV4", "Turbo15", "Turbo16", "Turbo17",
                      "Turbo18", "pATurbo15", "ApTurbo15", "pATurbo16",
                      "ApTurbo16", "Electron15", "Electron16", "Electron17",
                      "Electron18")
    if StripVer == '23':
        msg = ("""
            Attention: Stripping 23 is now replaced by Stripping 23_MCTuneV1. Please specify 23_MCTuneV1 as sample version."""
               )
        raise TypeError(msg)
    if StripVer == '23_MCTuneV1' or StripVer == '26' or StripVer == 'Turbo15' or StripVer == 'Turbo16' or StripVer == 'Turbo17' or StripVer == 'Turbo18' or StripVer == 'pATurbo15' or StripVer == 'ApTurbo15' or StripVer == 'pATurbo16' or StripVer == 'ApTurbo16' or StripVer == 'Electron15' or StripVer == 'Electron16' or StripVer == 'Electron17' or StripVer == 'Electron18':
        print(
            "***************************************************************************************************************************"
        )
        print(
            "INFO : ProbNN versions MC12TuneV2 and MC12TuneV3 should not be used for Run 2 data - advise to use MC12TuneV4 or MC15TuneV1"
        )
        print(
            "***************************************************************************************************************************"
        )
    if StripVer not in ValidStripVers:
        msg = ("Invalid stripping version '{0}'. "
               "Allowed versions are {1}").format(StripVer,
                                                  str(ValidStripVers))
        raise TypeError(msg)


def CheckStripVerPartNameMagPol(StripVer, PartName, MagPol):
    Valid_pA_Ap_PartNames = ["K", "P", "Pi", "Mu"]
    Valid_Electron_PartNames = ["e_B_Jpsi"]
    if 'pA' in StripVer or 'Ap' in StripVer:
        if 'Up' in MagPol:
            msg = (
                "MagUp data does not exist in the pA/Ap dataset. Please select MagDown and run again."
            )
            raise TypeError(msg)
        if PartName not in Valid_pA_Ap_PartNames:
            msg = (
                "Invalid choice of sample for pA/Ap data. Please use one of the following: K, Pi, Mu, P"
            )
            raise TypeError(msg)
    if 'Electron' in StripVer:
        if PartName not in Valid_Electron_PartNames:
            msg = (
                "Invalid choice of sample for 2015-2018 electron data. Please use the following: e_B_Jpsi"
            )
            raise TypeError(msg)
    #Check if the user has selected Turbo15 and Turbo16 along with e_B_Jpsi. Warn them to use Electron15 or Electron16 instead.
    if 'Turbo15' in StripVer or 'Turbo16' in StripVer or 'Turbo17' in StripVer or 'Turbo18' in StripVer:
        if PartName in Valid_Electron_PartNames:
            msg = (
                "The electron samples in Turbo15/16/17/18 should NOT be used! Please switch to using Electron15/16/17/18"
            )
            raise TypeError(msg)


def CheckVarName(VarName):
    ValidVarNames = GetVarNames()
    if VarName not in ValidVarNames:
        msg = ("Invalid binning variable '{0}'. "
               "Allowed variables are {1}").format(VarName, str(ValidVarNames))
        raise TypeError(msg)


def CheckIsTurbo(StripVer):
    ValidTurboNames = ("Turbo15", "Turbo16", "Turbo17", "Turbo18", "pATurbo15",
                       "ApTurbo15", "pATurbo16", "ApTurbo16")
    if StripVer not in ValidTurboNames:
        return False
    return True


#def GetTurboPath(StripVer,MagPol):
#    CheckStripVer(StripVer)
#    CheckMagPol(MagPol)
#    if StripVer=='Turbo15':
#        if MagPol=='MagDown':
#           return "/LHCb/Collision15/Beam6500GeV-VeloClosed-MagDown/Real Data/Reco15a/Turbo02/PIDCalibTuples2r1/95100000/PIDCALIB.ROOT"
#        if MagPol=='MagUp':
#           return "/LHCb/Collision15/Beam6500GeV-VeloClosed-MagUp/Real Data/Reco15a/Turbo02/PIDCalibTuples2r1/95100000/PIDCALIB.ROOT"
