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
#Added BKK functionality for retrieving Run 2 WGP nTuples
from DIRAC.Core.Base import Script
Script.initialize()
from DIRAC.Resources.Storage.StorageElement import StorageElement

import DIRAC
from LHCbDIRAC.BookkeepingSystem.Client.BookkeepingClient import BookkeepingClient
from DIRAC.DataManagementSystem.Client.DataManager import DataManager

#Read in the magnet polarity and the data version, along with the user path for writing the output to file
import sys
StripVer = sys.argv[1]
MagPolarity = sys.argv[2]
PartName = sys.argv[3]
path = sys.argv[4]

rm = DataManager()
bkClient = BookkeepingClient()
bkDict = {}
#args = Script.getPositionalArgs()
#print args

bkDict['ConfigName'] = 'LHCb'
bkDict['EventTypeId'] = '95100000'
bkDict['FileType'] = 'PIDCALIB.ROOT'
bkDict['DataTakingConditions'] = 'Beam6500GeV-VeloClosed-' + MagPolarity
#Proton-Proton collision data
if StripVer == "Turbo15" or StripVer == "Turbo16" or StripVer == "Turbo17" or StripVer == "Turbo18":
    year = StripVer.replace('Turbo', '')
    bkDict['ConfigVersion'] = 'Collision' + year
elif StripVer == "pATurbo15" or StripVer == "pATurbo16":
    year = StripVer.replace('pATurbo', '')
    bkDict['ConfigVersion'] = 'Protonion' + year
elif StripVer == "ApTurbo15" or StripVer == "ApTurbo16":
    year = StripVer.replace('ApTurbo', '')
    bkDict['ConfigVersion'] = 'Ionproton' + year

reco = ''
turbo = ''
version = ''
merge = ''
if StripVer == "Turbo15" or StripVer == "Turbo16" or StripVer == "Turbo17" or StripVer == "Turbo18":
    ############
    ### 2015 ###
    ############
    if year == '15':
        reco = '15a'
        turbo = '02'
        version = '5r1'
        merge = '05'
    ############
    ### 2016 ###
    ############
    elif year == '16':
        reco = '16'
        turbo = '02a'
        version = '9r3'
        #Jpsinopt files reproduced in dedicated production
        if (PartName == "Mu_nopt"):
            version = '9r3p1'
        merge = '05'
    ############
    ### 2017 ###
    ############
    elif year == '17':
        reco = '17'
        turbo = '04'
        version = '9r1'
        #Jpsinopt files reproduced in dedicated production
        if (PartName == "Mu_nopt"):
            version = '9r2p1'
        merge = '05'
    ############
    ### 2018 ###
    ############
    elif year == '18':
        reco = '18'
        turbo = '05'
        version = '9r2'
        #Jpsinopt files reproduced in dedicated production
        if (PartName == "Mu_nopt"):
            version = '9r2p1'
        merge = '05'

elif StripVer == "pATurbo15" or StripVer == "pATurbo16" or StripVer == "ApTurbo15" or StripVer == "ApTurbo16":
    ############
    ### 2015 ###
    ############
    if year == '15':
        reco = '15pLead'
        turbo = '03pLead'
        version = '5r0'
        merge = '01'
    ############
    ### 2016 ###
    ############
    elif year == '16':
        reco = '16pLead'
        turbo = '03pLead'
        version = '5r1'
        merge = '05'

bkDict[
    'ProcessingPass'] = '/Real Data/Reco' + reco + '/Turbo' + turbo + '/PIDCalibTuples' + version + '/PIDMerge' + merge

file = bkClient.getFiles(bkDict)
#file = bkClient.getFilesWithGivenDataSets(bkDict)
#file = bkClient.getVisibleFilesWithMetadata(bkDict)
#print file
files = file['Value']

print("There are " + str(len(files)) + " WGP nTuples in this dataset")
#print files
#print "============================================="
newFileList = []
for fileName in files:
    newFileList += [
        fileName.replace(
            '/lhcb/LHCb/',
            "root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/LHCb/")
    ]
    #print " ".join(newFileList)

#Write the file list to a json file for reading back in
import json

with open(path + "/files.txt", 'wt') as fp:
    json.dump(newFileList, fp)
