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

import os


#class _OCDSvars(object):
class _OCDSvars:
    _fileList = None
    _envVar = 0  #0 means unchecked, None means non-existent


def _warningMessage():
    print("")
    print(30 * '*' + " ! ! ! WARNING ! ! ! " + 29 * '*')
    print(
        "! You are using NON DEFAULT SAMPLES as the OVERRIDECALIBDATASTORE was set"
    )
    print("! Please make sure you understand why this warning appears.")
    print(80 * '*')
    print("")
    input("Press ENTER if you wish to continue...")


def _defineFileList():
    #  _warningMessage()
    _OCDSvars._fileList = []
    f = file(_getEnvVar())
    if not f:
        raise Exception("Error: File " + f + " not found.")

    for entry in f:
        _OCDSvars._fileList += [entry.replace('\n', '').replace('\r', '')]


def _getFileList():
    if _OCDSvars._fileList == None:
        _defineFileList()
    return _OCDSvars._fileList


def _getEnvVar():
    if _OCDSvars._envVar == 0:
        _OCDSvars._envVar = os.getenv("OVERRIDECALIBDATASTORE")
    return _OCDSvars._envVar


def _isEnvVarDefined():
    return _getEnvVar() != None and _getEnvVar() != ""


def GetDictFiles(runMin, runMax, maxFiles, verbose):
    if not _isEnvVarDefined(): return []

    if maxFiles == None or maxFiles < 0:
        maxFiles = 100000

    maxIndex = len(_getFileList())
    if maxIndex == None or int(maxIndex) > int(maxFiles):
        maxIndex = maxFiles

    return _getFileList()[:int(maxIndex)]


def GetDictIndex(runMin, runMax, maxFiles, verbose):
    if not _isEnvVarDefined(): return None

    if maxFiles == None or maxFiles < 0:
        maxFiles = 100000

    maxIndex = len(_getFileList())
    if maxIndex == None or maxIndex > maxFiles:
        maxIndex = maxFiles

    print("Min Index: 0;  Max Index: " + str(maxIndex - 1))

    ret = {'minIndex': 0, 'maxIndex': maxIndex - 1}
    return ret


def GetFileName(index):
    if not _isEnvVarDefined(): return None
    print("Required index: " + str(index))
    return _getFileList()[index]


if __name__ == '__main__':
    print("GetDictIndex ( 0, 100, 10, True )")
    print(GetDictIndex(0, 100, 10, True))
    print("-----------")
    print("GetFileName ( 1 )")
    print(GetFileName(1))
