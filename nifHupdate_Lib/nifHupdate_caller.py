#!/usr/bin/python
__author__="Claudia Chu"
__date__ ="2/18/19"

from nifHupdate_lib import parseConfig, createShFile, parseDate, launch, minimap_stages, throwError, testPrintFile

from os.path import abspath, join, isfile, isdir

from optparse import OptionParser

from os import curdir, system

from sys import argv, stderr, stdin

from sys import stdout

import time

import subprocess





#==============================================================
def real_main():
    usage  = "usage: %prog [options]"

    parser = OptionParser(usage)

    restartJob_def = 'minimap'
    parser.add_option( "-s", \
                       "--restartJob", \
                       type    = 'str', \
                       help    = "Start at this point in the pipeline.", \
                       default = restartJob_def )

    parser.set_description( "This pipeline requires a config file and name of a log file specified as arguments" )


    # Parsing the arguments
    (options, args) = parser.parse_args()

    # Moving to current directory
    basePath = abspath( curdir )

    nextStep = options.restartJob


    # Check arguments
    try:
        configFileName = argv[1]
    except:
        print("Please provide the path to config file")
        assert False
    #####

    try:
        logFile = argv[2]

    except:
        print("ERROR: Please provide name of log file")
        assert False
    #####

    # Setting of Configuration
    configFile = "%s/%s" % (basePath, configFileName)
    logFileFh = open(logFile, "a")
    configDict = parseConfig(configFile, basePath, logFileFh)

    if nextStep not in minimap_stages:
        print("'%s' is not a valid stage for minimap pipeline" % nextStep)
        assert False
    #####


    # Command list
    cmdList = []

    # Check if directory already exists. Prevent accidental overriding
    confirmDir = True
    if isdir(configDict["PREFIX"]):
        confirm = input("%s is already a directory. Are you sure you want to continue with this file name? (y/n) " % configDict["PREFIX"])
        confirmDir = confirm == "y"
        while (confirm != "y" and confirm != "n"):
            confirm = input("Please type either y for yes or n for no to continue using %s: " % configDict["PREFIX"]).strip()
            confirmDir = "y" == confirm
        #####
        if (not confirmDir):
            print("Cancelled")
            assert False
        #####
    else:
        if (confirmDir):
            newFolder = "mkdir %s" % (configDict["PREFIX"])
            cmdList.append(newFolder)
        #####
    #####


    cdFolder  = "cd %s/%s" % (basePath, configDict["PREFIX"])
    stageEcho = "echo Launching stage %s" % nextStep
    stageWrite = "echo Launching stage %s >> %s" % (nextStep, logFile)

    launchCmd = "python3 %s/nifHupdate_Lib/nifHupdate_launch.py %s %s %s %s/%s" % (basePath, configFile, nextStep, basePath, basePath, logFile)
    cdOut = "cd %s" % basePath # go back out
    cleanUpCmd = "mv %s_*.sh ./%s" % (configDict["PREFIX"], configDict["PREFIX"]) # Move stuff into directory

    cmdList.append("echo =========== RESTART ============ >> %s" % logFile)
    cmdList.append(stageEcho)
    cmdList.append(stageWrite)
    cmdList.append(cdFolder)
    cmdList.append(launchCmd)
    cmdList.append(cdOut)
    cmdList.append(cleanUpCmd)

    shFileName = createShFile(cmdList, basePath, configDict['PREFIX'], 'launch')
    # Stop if shFileName fails
    if ( not isfile(shFileName)):
        throwError("%s not available." % ( shFileName ))

    launch(shFileName)

#==============================================================
if ( __name__ == '__main__' ):
    real_main()

