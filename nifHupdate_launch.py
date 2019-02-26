#!/usr/bin/python
__author__="Claudia Chu"
__date__ ="2/24/19"

from nifHupdate_lib import parseConfig, createShFile, parseDate, launch, esearchCmds, fastaCmds

from os.path import abspath, join, isfile

from optparse import OptionParser

from os import curdir, system

from sys import stdout

import time

import subprocess

from sys import argv, stderr

# ======================================================

def real_main():
    # Read the command line information
    configFile  = argv[1]
    stage       = argv[2]
    basePath    = abspath( curdir )

    # Read the config file
    configDict = parseConfig( configFile )

    cmdList = []

    # parse years
    configDict = parseConfig(configFile)

    startDate = parseDate(configDict["START"])
    endDate = parseDate(configDict["END"])

    print(startDate)

    # for year in range(startDate, endDate):
    #     if (year == startDate):
    #         startDate = configDict["START"]
    #     if (year == endDate):
    #         startDate = configDict["END"]


# ================= STAGE 1 ===================== #
# esearch: getting the intial records and accession number

    if (stage == 'esearch'):
        print("esearch stage")
        for year in range (startDate, endDate):
            cmdList.append(esearchCmds(configDict, year))

        # Next stage
        nextStage = 'fasta'
        nextCmd = "python3 nifHupdate_launch.py %s %s" % (configFile, nextStage)
        cmdList.append(nextCmd)

        shFileName = createShFile(cmdList, basePath, configDict["PREFIX"], stage)
        launch(shFileName)





# ================= STAGE 2 ===================== #
    elif (stage == 'fasta'):
        print("fasta stage")
        for year in range (startDate, endDate):
            cmdList.append(fastaCmds(configDict, year, 'nifH'))
            cmdList.append(fastaCmds(configDict, year, 'genome'))

        # Next stage
        nextStage = 'lenCluster'
        nextCmd = "python3 nifHupdate_launch.py %s %s" % (configFile, nextStage)
        cmdList.append(nextCmd)
        shFileName = createShFile(cmdList, basePath, configDict["PREFIX"], stage)
        launch(shFileName)

    else:
        print("Done!")


#==============================================================
if ( __name__ == '__main__' ):
    real_main()