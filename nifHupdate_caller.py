#!/usr/bin/python
__author__="Claudia Chu"
__date__ ="2/18/19"

from nifHupdate_lib import parseConfig, createShFile, parseDate, launch

from os.path import abspath, join, isfile

from optparse import OptionParser

from os import curdir, system

from sys import stdout

import time

import subprocess

#==============================================================
def setParser():
    usage  = "usage: %prog [options]"

    parser = OptionParser(usage)

    restartJob_def = 'megahit'
    parser.add_option( "-s", \
                       "--restartJob", \
                       type    = 'str', \
                       help    = "Start at this point in the pipeline.", \
                       default = restartJob_def )

    parser.set_description( "This pipeline requires a config file called 'nifHupdate.config' in the current working directory. It must include" )


    # Parsing the arguments
    (options, args) = parser.parse_args()



#==============================================================
def real_main():
    setParser()
     # Moving to current directory
    basePath = abspath( curdir )

    # Setting the config file name
    configFile = "%s/nifHupdate_config.txt" % basePath

    configDict = parseConfig( configFile )

    launch_cmd = "python3 nifHupdate_launch.py %s %s" % (configFile, 'esearch')


    shFileName = createShFile([launch_cmd], basePath, configDict['PREFIX'], 'launch')
    if ( not isfile(shFileName)):
        throwError("%s not available" % ( shFileName ) )
    launch(shFileName)

#==============================================================
if ( __name__ == '__main__' ):
    real_main()

