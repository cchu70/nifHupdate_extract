#!/usr/bin/python
__author__="Claudia Chu"
__date__ ="2/18/19"

from nifHupdate_lib import parseConfig, createShFile, parseDate, launch, stages, throwError

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
        throwError("%s is not a valid stage." % nextStep)
    #####


    configFileName = argv[1]

    # Setting the config file name
    configFile = "%s/%s" % (basePath, configFileName)

    configDict = parseConfig(configFile, basePath)
    if isdir(configDict["PREFIX"]):
        print("%s had already been used. Are you sure you want to continue with this file name?", configDict["PREFIX"])
        confirm = stdin
        if (confirm == "n"):
            print("Cancelled")
            assert False
        #####
    #####



    cmdList = []
    newFolder = "mkdir %s" % (configDict["PREFIX"])
    cdFolder  = "cd %s" % (configDict["PREFIX"])
    launchCmd = "python3 %s/nifHupdate_launch.py %s %s %s" % (basePath, configFile, nextStep, basePath + configDict["PREFIX"])


    cmdList.append(newFolder)
    cmdList.append(cdFolder)
    cmdList.append(launchCmd)

    shFileName = createShFile(cmdList, basePath, configDict['PREFIX'], 'launch')
    if ( not isfile(shFileName)):
        throwError("%s not available." % ( shFileName ) )
    launch(shFileName)

#==============================================================
if ( __name__ == '__main__' ):
    real_main()

