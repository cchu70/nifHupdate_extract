#!/usr/bin/python
__author__="Claudia Chu"
__date__ ="2/24/19"

from nifHupdate_lib import parseConfig, createShFile, parseDate, launch, esearchCmds, fastaCmds, blastnCmds, bestAlignment, throwError, verifyDb

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
    basePath    = argv[3]

    # Read the config file
    configDict = parseConfig(configFile, basePath)

    cmdList = []

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
        nextCmd = "python3 %s/nifHupdate_launch.py %s %s %s\n" % (basePath, configFile, nextStage, basePath)
        cmdList.append(nextCmd)

        shFileName = createShFile(cmdList, basePath, configDict["PREFIX"], stage)
        launch(shFileName)


# ================= STAGE 2 ===================== #
    elif (stage == 'fasta'):
        print("fasta stage")
        for year in range (startDate, endDate):
            for sortterm in configDict["SORTTERMS"]:
                fastaFileName = "%s_%s.%s.fasta" % (configDict["PREFIX"], year, sortterm)
                cmdList.append(fastaCmds(configDict, year, sortterm, fastaFileName))
            #####
        #####

        # Next stage
        nextStage = 'set_db'
        nextCmd = "python3 %s/nifHupdate_launch.py %s %s %s\n" % (basePath, configFile, nextStage, basePath)
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
        nextCmd = "python3 %s/nifHupdate_launch.py %s %s %s\n" % (basePath, configFile, nextStage, basePath)
        cmdList.append(nextCmd)

        # n = subprocess.Popen()


# ================= STAGE 4 ===================== #
# Aligning each of the query sequences to the subject sequences from local database.

    elif (stage == 'blastn'):
        print("blastn")

        fofnFileName = "%s.blastnFiles.fofn" % (configDict["PREFIX"])
        fh = open(fofnFileName, "w")

        for year in range (startDate, endDate):
            for sortterm in configDict["SORTTERMS"]:
                cmdList.append(blastnCmds(configDict, year, sortterm, fh))

        fh.close()

        nextStage = 'filter_best_alignments'
        nextCmd = "python3 %s/nifHupdate_launch.py %s %s %s\n" % (basePath, configFile, nextStage, basePath)

        cmdList.append(nextCmd)
        shFileName = createShFile(cmdList, basePath, configDict["PREFIX"], stage)
        launch(shFileName)

# ================= STAGE 5 ===================== #
# Parsing through blastn tables to find best alignment for each sequence

    elif (stage == 'filter_best_alignments'):


        fofnFileName = "%s.blastnFiles.fofn" % (configDict["PREFIX"])
        if (not isfile(fofnFileName)) :
            throwError("%s is not available. Either re-run blastn stage, or cat all your blastn files into this file name." % fofnFileName)
        #####

        for sortterm in configDict["SORTTERMS"]:
            fileName = "%s.%s.blastn.txt" % (configDict["PREFIX"], sortterm)
            fh = open(fileName, "w")
            for blastnFile in open(fofnFileName, "r"):
                if sortterm in blastnFile:
                    bestAlignment(blastnFile.strip(), fh)
                #####
            #####
            fh.close()
        #####

        nextStage = 'trim_seq'
        nextCmd = "python3 %s/nifHupdate_launch.py %s %s %s\n" % (basePath, configFile, nextStage, basePath)

        cmdList.append(nextCmd)
        shFileName = createShFile(cmdList, basePath, configDict["PREFIX"], stage)
        launch(shFileName)

# ================= STAGE 6 ===================== #

    elif (stage == 'trim_seq'):

        nextStage = 'end'
        nextCmd = "python3 %s/nifHupdate_launch.py %s %s %s\n" % (basePath, configFile, nextStage, basePath)

        # parse the blastn files into a map
        blastnFofn = "%s.blastn.fofn" % (configDict["PREFIX"])
        blastnMap = mapBlast(blastnFofn)
        # parse the esearch files into a map
        esearchFofn = "%s.esearch.fofn" % (configDict["PREFIX"])
        esearchMap = mapEsearch(esearchFofn)


        # for each year and each fasta file (fa_list.txt)
        # 1) get the accession number
        # 2) get associated blastn data: accession number, cluster, type?, cds/genome
        # 3) make new record name with updated data
        # 4) trim query sequence (qstart qend)
        # 5) print to a new file

        # >AF484654;cluster_I;Mesorhizobium_sp._LMG_11892

        for fastaFile in open("fa_list.txt", "r"):


        cmdList.append(nextCmd)
        shFileName = createShFile(cmdList, basePath, configDict["PREFIX"], stage)
        launch(shFileName)

# ================= STAGE 7 ===================== #
    #elif (stage == 'cluster'):
#==============================================================
if ( __name__ == '__main__' ):
    real_main()