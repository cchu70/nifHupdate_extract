#!/usr/bin/python
__author__="Claudia Chu"
__date__ ="2/24/19"

from nifHupdate_lib import parseConfig, createShFile, parseDate, launch, esearchCmds, fastaCmds, blastnCmds, bestAlignment, throwError, verifyDb, test, testPrintFile, wait, fasta

from os.path import abspath, join, isfile

from optparse import OptionParser

from os import curdir, system, stat

from sys import stdout

from Bio import SeqIO

import time

import subprocess

from sys import argv, stderr

import threading

# ======================================================

def real_main():

    # TESTING
    # test("In Launch File")
    # ######


    # Read the command line information
    configFile  = argv[1]
    stage       = argv[2]
    basePath    = argv[3]

    # TESTING
    # test("ConfigFile %s" % configFile.split("/")[-1], "stage %s" % stage, "basePath %s" % basePath)
    ######

    # Read the config file
    configDict = parseConfig(configFile, basePath)

    # Get numerical year
    startDate = parseDate(configDict["START"])
    endDate = parseDate(configDict["END"])

    # for year in range(startDate, endDate):
    #     if (year == startDate):
    #         startDate = configDict["START"]
    #     if (year == endDate):
    #         startDate = configDict["END"]


    # FREQUENTLY USED CONSTANTS AND VARIABLES
    CMDLIST = []
    PREFIX = configDict["PREFIX"]


    # # Move into current running directory
    # cdFolder  = "cd %s" % (configDict["PREFIX"])
    # CMDLIST.append(cdFolder)

    # TESTING
    # test("PREFIX: %s" % PREFIX, "STARTDATE: %s" % startDate, "ENDDATE: %s" % endDate)
    ######


# ================= STAGE 1 ===================== #
# esearch: getting the intial records and accession number

    if (stage == 'esearch'):

        # TESTING
        # test("In Esearch!")
        ######

        # Track esearch text files
        esearchFofn = "%s.esearch.fofn" % PREFIX
        fh = open(esearchFofn, "w")



        for year in range (startDate, endDate):
            outputfile = "%s_%s.esearch.txt\n" % (PREFIX, year)

            # Replace with literal year

            if (year == startDate):
                year = configDict["START"]
            elif (year == endDate):
                year = configDict["END"]
            #####

            CMDLIST.append("echo esearch for year %s..." % year)
            CMDLIST.append(esearchCmds(configDict, year, outputfile))

            fh.write(outputfile)

        fh.close()
        # Next stage
        nextStage = 'fasta'
        echoCmd = "echo Launching next stage %s" % nextStage
        nextCmd = "python3 %s/nifHupdate_Lib/nifHupdate_launch.py %s %s %s\n" % (basePath, configFile, nextStage, basePath)
        CMDLIST.append(echoCmd)
        CMDLIST.append(nextCmd)

        shFileName = createShFile(CMDLIST, basePath, PREFIX, stage)


        # # TESTING
        # testPrintFile(shFileName)
        # ######

        launch(shFileName)


# ================= STAGE 2 ===================== #
    elif (stage == 'fasta'):

        # # TESTING
        # test("In fasta!")
        # ######

        fastaFofn = "%s.fasta.fofn" % PREFIX
        fh = open(fastaFofn, "w")
        esearchFofn = "%s/%s/%s.esearch.fofn" % (basePath, PREFIX, PREFIX)
        if (not isfile(esearchFofn)):
            throwError("esearch stage failed: %s not found" % esearchFofn)
        else:


            threads = []
            for esearchFile in open(esearchFofn, "r"):
                if (not isfile(esearchFile.strip())):
                    throwError("esearch stage failed: %s not found" % file)
                elif (not stat(esearchFile.strip())):
                    throwError("esearch stage failed: %s empty" % file)
                #####

                # Make the fasta commands
                prefix = esearchFile.split(".")[0]
                for sortterm in configDict["SORTTERMS"]:
                    fastaFileName = "%s.%s.fasta" % (prefix, sortterm)
                    fh.write("%s\n" % fastaFileName) # record fasta file names

                    print("Retrieving %s" % fastaFileName)

                    # # OLD
                    # CMDLIST.append("echo Retrieving %s ..." % fastaFileName)
                    # CMDLIST.append(fastaCmds(esearchFile.strip(), sortterm, basePath, fastaFileName))

                    t = threading.Thread(target=fasta, args= (esearchFile.strip(), sortterm, fastaFileName))
                    threads.append(t)
                    t.start()
                    t.join()
                    #fasta(esearchFile.strip(), sortterm, fastaFileName)
                #####
            #####
        #####
        # # wait for threads to finish
        # while (not t.isAlive() for t in threads):
        #     time.sleep(5);
        #     print("still alive")

        fh.close()


        # Next stage
        nextStage = 'set_db'
        nextCmd = "python3 %s/nifHupdate_Lib/nifHupdate_launch.py %s %s %s\n" % (basePath, configFile, nextStage, basePath)
        CMDLIST.append(nextCmd)
        shFileName = createShFile(CMDLIST, basePath, PREFIX, stage)

        # # TESTING
        # testPrintFile(shFileName)
        # ######
        launch(shFileName)

# ================= STAGE 3 ===================== #
# Check if all the files exist for a local database for stand alone alignment

    elif (stage == 'set_db'):
        # TESTING
        test("In set_db!")
        ######

        # check database already exists
        if (not verifyDb(configDict["DBFILE"])):
            # move to next stage
            makeDbCmd = "makeblastdb -in %s/%s -parse_seqids -dbtype nucl -out %s/%s/%s" % (basePath, configDict["DBFILE"], basePath, PREFIX, configDict["DBNAME"])
            CMDLIST.append(makeDbCmd)
        #####

        nextStage = 'blastn'
        nextCmd = "python3 %s/nifHupdate_Lib/nifHupdate_launch.py %s %s %s" % (basePath, configFile, nextStage, basePath)
        CMDLIST.append(nextCmd)

        shFileName = createShFile(CMDLIST, basePath, PREFIX, stage)
        # # TESTING
        # testPrintFile(shFileName)
        # ######
        launch(shFileName)


# ================= STAGE 4 ===================== #
# Aligning each of the query sequences to the subject sequences from local database.

    elif (stage == 'blastn'):
        # # TESTING
        # test("In blastn!")
        # ######

        fastaFofn = "%s.fasta.fofn" % PREFIX
        if (not isfile(fastaFofn)):
            throwError("fasta stage failed: %s not found" % esearchFofn)
        else:

            fofnFileName = "%s.blastnFiles.fofn" % (PREFIX)
            fh = open(fofnFileName, "w")

            for fastaFile in open(fastaFofn, "r"):
                prefix, source, end = fastaFile.strip().split(".")
                outputFile = "%s.%s.blastn.txt" % (prefix, source)
                fh.write(outputFile + "\n") # write to blast fofn file
                CMDLIST.append("echo Blasting %s ..." % fastaFile.strip())
                CMDLIST.append(blastnCmds(configDict["DBNAME"], fastaFile.strip(), outputFile))

            fh.close()
        #####

        nextStage = 'filter_best_alignments'
        nextCmd = "python3 %s/nifHupdate_Lib/nifHupdate_launch.py %s %s %s\n" % (basePath, configFile, nextStage, basePath)

        CMDLIST.append(nextCmd)
        shFileName = createShFile(CMDLIST, basePath, PREFIX, stage)

        # TESTING
        testPrintFile(shFileName)
        ######

        launch(shFileName)

# ================= STAGE 5 ===================== #
# Parsing through blastn tables to find best alignment for each sequence

    elif (stage == 'filter_best_alignments'):
        # TESTING
        test("In filter_best_alignments!")
        ######

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
        nextCmd = "python3 %s/nifHupdate_Lib/nifHupdate_launch.py %s %s %s\n" % (basePath, configFile, nextStage, basePath)

        CMDLIST.append(nextCmd)
        shFileName = createShFile(CMDLIST, basePath, configDict["PREFIX"], stage)
        launch(shFileName)

# ================= STAGE 6 ===================== #

    elif (stage == 'trim_seq'):

        assert False
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

        # for fastaFile in open("fa_list.txt", "r"):
        fastaFofn = "%s.fasta.fofn" % PREFIX
        #for fastaFile in open(fastaFofn, "r"):



        CMDLIST.append(nextCmd)
        shFileName = createShFile(CMDLIST, basePath, configDict["PREFIX"], stage)
        launch(shFileName)

# ================= STAGE 7 ===================== #
    #elif (stage == 'cluster'):
# ================= STAGE 8 ===================== #
    #elif (stage == 'duplicate'):
    # cd-hit-dup -i fasta -o output


    # Reading stderr


    return 0
#==============================================================
if ( __name__ == '__main__' ):
    real_main()