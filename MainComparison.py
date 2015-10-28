import os,csv,shutil,sys
import string
import sys
import csv
import os
import shutil
import time, argparse,logging,subprocess
import logging.handlers,re
import ConfigParser
from itertools import izip


fileType = "txt"

def logging_setup(logger,logpath,jobName):
    if not os.path.exists(logpath) or not os.path.isdir(logpath):
        logpath = "./"
    log_path = os.path.join(logpath, "./LOGGING.txt")
    #log_warn_path = os.path.join(logpath,"./" + jobName + "_WARNINGS.txt")
    #log_error_path = os.path.join(logpath, "./" + jobName + "_ERRORS.txt")
    formatter = logging.Formatter('[%(asctime)s] %(levelname)-8s: %(message)s', datefmt='%m-%d %H:%M')
    infofh = logging.handlers.TimedRotatingFileHandler(os.path.abspath(log_path))
    infofh.setLevel(logging.INFO)
    infofh.setFormatter(formatter)
    warningfh = logging.handlers.TimedRotatingFileHandler(os.path.abspath(log_warn_path))
    warningfh.setLevel(logging.WARN)
    warningfh.setFormatter(formatter)
    errorfh = logging.handlers.TimedRotatingFileHandler(os.path.abspath(log_error_path))
    errorfh.setLevel(logging.ERROR)
    errorfh.setFormatter(formatter)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    logger.setLevel(logging.INFO)
    logger.addHandler(ch)
    logger.addHandler(infofh)
    logger.addHandler(warningfh)
    logger.addHandler(errorfh)
    return



def parseData(fileName):
    print "Parsing data for " + fileName
    resultDict = dict()
    with open(fileName,mode = 'r') as inF:
        geneName = ""
        currString = ""
        firstGene = True
        for line in inF:
            if ">" in line:
                if firstGene:
                    #print "at first gene"
                    firstGene = False
                else:
                   # print "assigning gene "
                    resultDict[geneName] = currString
                    geneName = ""
                    currString = ""
                    
                geneLineList = line.split()
                geneLine = geneLineList[2] #sequence of comma seperated names is third item in list
                geneSplitList = geneLine.split(',') #gene name is the first string in the sequence of comma seperated variables
                geneName = str(geneSplitList[0])
            else:
                if line == "\n":
                    pass
                elif "\n" in line:
                    currString += line[:len(line)-1] #gets rid of \n charecter
                else:
                    currString += line[:len(line)]
        resultDict[geneName] = currString
    inF.close()

    for key in resultDict:
        #if key == "YBR160W":
        #print key,resultDict[key]
        pass
    #print resultDict
    return resultDict

def parseMain(fileName):
    print "Parsing data for " + fileName
    resultDict = dict()
    firstGene = True
    with open(fileName,mode = 'r') as inF:
        geneName = ""
        currString = ""
        for line in inF:
            #means its a gene name
            #print len(line)
            #print repr(line)
            if ">" in line:
                """
                if firstGene:
                    #If its the first gene, don't write to dict yet
                    firstGene = False
                else:
                    print "gosu"
                    resultDict[geneName] = currString
                    currString,geneName = "",""
                """
                geneName,currString = "",""
                geneLineList = line.split()
                geneLine = geneLineList[0]
                geneName = geneLine[1:] #gets rid of the '<' at the beginning of the name
                resultDict[geneName] = ""
                #print "found a gene: " + geneName
                #print resultDict[geneName]
            #else its a sequence so add to current gene name
            elif len(line) == 61:
                #print "found a line of genes: " + line
                if "\n" in line:
                    currString += line[:len(line)-1]
                else:
                    currString += line[:len(line)]
            else:
                assert(">" not in line)
                if "\n" in line:
                    currString += line[:len(line)-1]
                else:
                    currString += line[:len(line)]
                resultDict[geneName] = currString
                currString,geneName = "",""
    inF.close()
    for key in resultDict:
        #if key == "YBR160W":
        #print key,resultDict[key]
        pass
    return resultDict

def compare(item,other):
    return item == other

#Check to see if there are any genes in the main strain that are not existent in the comparitive strain
def compareData(mainStrain,secondDict,thirdDict = None):
    firstDictKeys = mainStrain.keys()
    print "comparing data"
    differentGeneList = list()
    sameGeneList = list()
    resultFile = os.path.join(os.getcwd(),"WineResults.csv")
    otherFile = resultFile
    if os.path.isfile(resultFile):
        count = 1
        while os.path.isfile(resultFile):
            resultFile = otherFile[:len(otherFile)-4] + "(" + str(count) + ").csv"
            count +=1
    with open(resultFile, 'wb') as csvfile:
        spamwriter = csv.writer(csvfile,delimiter = ',')
        spamwriter.writerow(["Gene Name","Main Strain Length", "Other Strain Length", "Number of differences", "List of differences"])

        for key in mainStrain:
            #print key
            splitKeyList = key.split()
            keyName = splitKeyList[0]
            if keyName in secondDict:
                #Common genes
                mainStrainSeq = mainStrain[key]
                otherStrainSeq = secondDict[key]
                #sameGeneList.append(key)
                #print compare(mainStrainSeq,otherStrainSeq)
                differenceList = [i for i, (a1,a2) in enumerate(izip(mainStrainSeq,otherStrainSeq)) if a1!=a2]
                #spamwriter.writerow([keyName, len(mainStrainSeq),len(otherStrainSeq),len(differenceList),differenceList])
                spamwriter.writerow([keyName,len(mainStrainSeq),len(otherStrainSeq),len(differenceList),differenceList])
                if len(differenceList) >= 1:
                    sameGeneList.append(keyName)

            elif keyName not in secondDict:
                differentGeneList.append(key)
       # print differentGeneList
        #print len(sameGeneList)
       # print len(differentGeneList),len(sameGeneList),len(mainStrain),len(secondDict)
    pass

def extractData(mainStrain = None,secondFile = None, thirdFile = None):
    #firstResultDict = dict() #dictionary where key is gene and entry is DNA sequence
    #secondResultDict = dict()
    #print "files one and two are " + mainStrain + " " + secondFile
    mainStrainDict = dict()
    secondResultDict = dict()
    thirdResultDict = dict()
    if (mainStrain != None and mainStrain.endswith(fileType)):
        mainStrainDict = parseMain(mainStrain)
    
    #for keys in firstResultDict:
    #    print keys, firstResultDict[keys]


    if (secondFile != None and secondFile.endswith(fileType)):
        secondResultDict = parseData(secondFile)

    #if (thirdFile != None and thirdFile.endswith(fileType)):
    #    thirdResultDict = parseData(thirdFile)

    print len(mainStrainDict), len(secondResultDict)#,len(thirdResultDict) 
    #maxKeys = max(len(ResultDict),len(secondResultDict))#,len(thirdResultDict))
    
    if len(mainStrainDict) == len(secondResultDict):
        print "Both strains have same number of genes."
    else:
        print "They have different number of genes."
    
    mainDict = dict()
    compareData(mainStrainDict,secondResultDict)

"""
    if maxKeys == len(firstResultDict):
        compareData(firstResultDict,secondResultDict)
    elif maxKeys == len(secondResultDict):
        compareData(secondResultDict,firstResultDict)
#    elif maxKeys == len(thirdResultDict):
#        print "nops"
    else:
        #all same length
        print "gpsu"
"""

if __name__ == '__main__':
    #log = logging.getlogger("WineTest_logger")
    #logging_setup(log)
    if len(sys.argv) == 1:
        print "Cannot proceed, not enough input files..."
    elif len(sys.argv) == 2:
        print "Please provide one more input file..."
    else:
        assert(len(sys.argv) > 2)
        mainStrainFile = sys.argv[1]
        for i in xrange(len(sys.argv) - 2):
            extractData(mainStrainFile,sys.argv[i+2])


