#!/usr/bin/python
__author__="Claudia Chu"
__date__ ="2/18/19"

from nifHupdate_lib import parseConfig, createShFile, parseDate, launch, stages, throwError, testPrintFile

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

    restartJob_def = 'esearch'
    parser.add_option( "-s", \
                       "--restartJob", \
                       type    = 'str', \
                       help    = "Start at this point in the pipeline.", \
                       default = restartJob_def )

    parser.set_description( "This pipeline requires a config file specified as an argument" )


    # Parsing the arguments
    (options, args) = parser.parse_args()

    # Parsing the arguments
    (options, args) = parser.parse_args()

    # Moving to current directory
    basePath = abspath( curdir )

    nextStep = options.restartJob

    # Check if stage is valid
    if nextStep not in stages:
        throwError("'%s' is not a valid stage." % nextStep)
    #####

    try:
        configFileName = argv[1]
    except:
        print("Please provide the path to config file")
        assert False
    #####

    try:
        logFile = argv[2]
    except:
        print("Please provide the path to log file")
        assert False
    #####

    # Setting the config file name
    configFile = "%s/%s" % (basePath, configFileName)



    # Command list
    cmdList = []

    configDict = parseConfig(configFile, basePath)

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
            //throwError("Cancelled", )
        #####
    else:
        if (confirmDir):
            newFolder = "mkdir %s" % (configDict["PREFIX"])
            cmdList.append(newFolder)
        #####
    #####


    # Continue adding commands
    cdFolder  = "cd %s" % (configDict["PREFIX"])
    launchCmd = "python3 %s/nifHupdate_Lib/nifHupdate_launch.py %s %s %s %s" % (basePath, configFile, nextStep, basePath, logFile)

    cmdList.append(cdFolder)
    cmdList.append(launchCmd)

    shFileName = createShFile(cmdList, basePath, configDict['PREFIX'], 'launch')
    # Stop if shFileName fails
    if ( not isfile(shFileName)):
        throwError("%s not available." % ( shFileName ) )

    # ####### TESTING ######
    # printFile(shFileName)
    # assert False
    # ######################
    launch(shFileName)

#==============================================================
if ( __name__ == '__main__' ):
    real_main()

