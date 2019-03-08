#!/usr/bin/python
__author__="Claudia Chu"
__date__ ="2/24/19"

from nifHupdate_lib import parseConfig, createShFile, parseDate, launch, esearchCmds, fastaCmds, blastnCmds, bestAlignment

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
        nextStage = 'blastn'
        nextCmd = "python3 nifHupdate_launch.py %s %s" % (configFile, nextStage)
        cmdList.append(nextCmd)
        shFileName = createShFile(cmdList, basePath, configDict["PREFIX"], stage)
        launch(shFileName)

# ================= STAGE 3 ===================== #
# Check if all the files exist for a local database for stand alone alignment

    elif (stage == 'set_db'):
        # check database already exists
        if (not verifyDb(configDict["DBNAME"])):
            # move to next stage
            makeDbCmd = "makeblastdb -in %s -parse_seqids -title 'TAXADIVA' -dbtype nucl" % (configDict["DBNAME"])
            cmdList.append(makeDbCmd)
        #####

        nextStage = 'blastn'
        nextCmd = "python3 nifHupdate_launch.py %s %s" % (configFile, nextStage)
        cmdList.append(nextCmd)

        # n = subprocess.Popen()


# ================= STAGE 4 ===================== #
# Aligning each of the query sequences to the subject sequences from local database.

    elif (stage == 'blastn'):
        print("blastn")

        fofnFileName = "%s.blastnFiles.fofn" % (configDict["PREFIX"])
        fh = open(fofnFileName, "w")

        for year in range (startDate, endDate):
            cmdList.append(blastnCmds(configDict, year, 'gene', fh))

        for year in range (startDate, endDate):
            cmdList.append(blastnCmds(configDict, year, 'genomes', fh))

        fh.close()

        nextStage = 'filter_alignments'
        nextCmd = """python3 nifHupdate_launch.py %s %s\n""" % (configFile, nextStage)

        cmdList.append(nextCmd)
        shFileName = createShFile(cmdList, basePath, configDict["PREFIX"], stage)
        launch(shFileName)

# ================= STAGE 5 ===================== #
# Parsing through blastn tables to find best alignment for each sequence

    elif (stage == 'filter_alignments'):


        fofnFileName = "%s.blastnFiles.fofn" % (configDict["PREFIX"])

        geneFiltered = "%s.%s.blastn.txt" % (configDict["PREFIX"], "gene")
        genomeFiltered = "%s.%s.blastn.txt" % (configDict["PREFIX"], "genome")
        # storing the correct alignment

        nh = open(geneFiltered, "w")
        gh = open(genomeFiltered, "w")

        for blastnFile in open(fofnFileName, "r"):
            if "gene" in blastnFile:
                bestAlignment(blastnFile.strip(), nh)
            else:
                bestAlignment(blastnFile.strip(), gh)
            #####
        #####

        nh.close()
        gh.close()

        nextStage = 'trim'
        nextCmd = """python3 nifHupdate_launch.py %s %s\n""" % (configFile, nextStage)

        cmdList.append(nextCmd)
        shFileName = createShFile(cmdList, basePath, configDict["PREFIX"], stage)
        launch(shFileName)

# ================= STAGE 6 ===================== #

    elif (stage == 'trim'):

        nextStage = 'end'
        nextCmd = """python3 nifHupdate_launch.py %s %s\n""" % (configFile, nextStage)

        cmdList.append(nextCmd)
        shFileName = createShFile(cmdList, basePath, configDict["PREFIX"], stage)
        launch(shFileName)

#==============================================================
if ( __name__ == '__main__' ):
    real_main()