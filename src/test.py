# coding: utf8
import DnaUtils
import FastqSequence
import logging
from logging.config import fileConfig
import os
import sys
import re

def tryToFusion(fastqSequence, reverseFileName, threshold, log, seqOutput):
    fastqSequence.applyDrasticThreshold(threshold)
    f = open(reverseFileName, "r")
    reverseLineHeader = fastqSequence.getReverseLineHeader()
    logging.debug("Search " + reverseLineHeader)
    
    lastHeaderLineNumber    = 0
    lastSeparatorLineNumber = 0
    parseFastq     = False
    rFastqSequence = None
    lineNumber     = 0
    
    for line in f:
        line = line[:-1]
        lineNumber += 1
        #logging.debug(str(lineNumber) +"\t****   "+line + "\t" + str(line[0] == "@") + "\t" + str(lineNumber == 1))
        
        if line[0] == "@" and (lastSeparatorLineNumber != lineNumber - 1 or lineNumber == 1):
            logging.debug(line + " est un header")
            if parseFastq:
                rFastqSequence.applyDrasticThreshold(threshold)
                
                s = FastqSequence.matchFastqSequence(fastqSequence, rFastqSequence)
                if s != None:
                    seqOutput.write(">" + s.getLineHeader() + "\n")
                    seqOutput.write(s.sequence + "\n")
                    logging.debug(s.getLineHeader() + "\t" + str(s.typeFusion))
                    if s.typeFusion == s.TYPE_FUSION_OK:
                        log.write("Merge " + fastqSequence.getLineHeader() + " with " + rFastqSequence.getLineHeader() + "\n")
                    elif s.typeFusion == s.TYPE_FUSION_NO_MATCH:
                        log.write("Forced merge " + fastqSequence.getLineHeader() + " with " + rFastqSequence.getLineHeader() + "\n")
                    elif s.typeFusion == s.TYPE_FUSION_PARTIAL:
                        log.write("Partial merge " + fastqSequence.getLineHeader() + " with " + rFastqSequence.getLineHeader() + "\n")
                    else:
                        log.write("Can't merge " + fastqSequence.getLineHeader() + " with " + rFastqSequence.getLineHeader() + "\n")
                else:
                    log.write("*Can't merge " + fastqSequence.getLineHeader() + " with " + rFastqSequence.getLineHeader() + "\n")
                parseFastq = False
                
            #logging.debug("-"+line)
            #logging.debug("-"+reverseLineHeader)
            #logging.debug(line + " =?= " + reverseLineHeader + " ? " + str(line == reverseLineHeader))
            if line == reverseLineHeader:
                parseFastq = True
                
            if parseFastq:
                rFastqSequence = FastqSequence.readFasqSequenceHeader(line)
                #logging.debug("Compare to " + rFastqSequence.getLineHeader())
                lastHeaderLineNumber = lineNumber
        elif parseFastq:
            if  lastHeaderLineNumber == lineNumber - 1:
                rFastqSequence.setSequence(line)
            elif line[0] == "+":
                lastSeparatorLineNumber = lineNumber
            elif lastSeparatorLineNumber == lineNumber - 1:
                rFastqSequence.setQuality(line)                
            
    f.close()

def testQualityThreshold(fastqFileName, threshold):
    f = open(fastqFileName, "r")
    
    lineNumber = 0
    lastHeaderLineNumber    = 0
    lastSeparatorLineNumber = 0
    sequences  = []
    fastqSequence = None
    
    for line in f:
        line = line[:-1]
        if line[0] == "@" and lastSeparatorLineNumber != lineNumber - 1:
            if fastqSequence != None:
                #tryToFusion(fastqSequence, reverseFileName, threshold, log, seqOutput)
                avLen = len(fastqSequence.sequence)
                fastqSequence.applyDrasticThreshold(threshold)
                logger.debug("Avant : " + str(avLen) + " / Apres : " + str(len(fastqSequence.sequence)))
                
            fastqSequence = FastqSequence.readFasqSequenceHeader(line)
            lastHeaderLineNumber = lineNumber
            
            logging.debug(line)
        elif  lastHeaderLineNumber == lineNumber - 1:
            fastqSequence.setSequence(line)
        elif line[0] == "+":
            lastSeparatorLineNumber = lineNumber
        elif lastSeparatorLineNumber == lineNumber - 1:
            fastqSequence.setQuality(line)
            
        lineNumber += 1
    #tryToFusion(fastqSequence, reverseFileName, threshold, log, seqOutput)
    f.close()

def readFastQFile(fileName, reverseFileName, threshold, seqOutputName, logFileName):
    f = open(fileName, "r")
    log = open(logFileName, "w")
    seqOutput = open(seqOutputName, "w")
    
    lineNumber = 0
    lastHeaderLineNumber    = 0
    lastSeparatorLineNumber = 0
    sequences  = []
    fastqSequence = None
    
    for line in f:
        line = line[:-1]
        if line[0] == "@" and lastSeparatorLineNumber != lineNumber - 1:
            if fastqSequence != None:
                tryToFusion(fastqSequence, reverseFileName, threshold, log, seqOutput)
                
            fastqSequence = FastqSequence.readFasqSequenceHeader(line)
            lastHeaderLineNumber = lineNumber
            
            logging.debug(line)
        elif  lastHeaderLineNumber == lineNumber - 1:
            fastqSequence.setSequence(line)
        elif line[0] == "+":
            lastSeparatorLineNumber = lineNumber
        elif lastSeparatorLineNumber == lineNumber - 1:
            fastqSequence.setQuality(line)
            
        lineNumber += 1
    tryToFusion(fastqSequence, reverseFileName, threshold, log, seqOutput)
    f.close()
    log.close()
    seqOutput.close()   
    
    #for s in sequences:
    #    s.applyDrasticThreshold(threshold)
    
    #f = FastqSequence.matchFastqSequence(sequences[0], sequences[1])
    #print(f.sequence)
    #print(sequences[0] == sequences[1])

#def searchSequence(fastqSequence, fileName):
    
def checkArgument():
    # print(DnaUtils.reverseSequence(DnaUtils.complementDnaSequence("GATATTTTAGATATTTTCCGAAGGTACCATTTAA")))
    # sys.exit(1)
    if len(sys.argv) -1 != 5:
        print(sys.argv[0] + " [R1.fastq] [R2.fastq] [drastic threshold] [seq ouput] [log ouputFile]")
        print("Exemple " + sys.argv[0] + " P1-A01_GACGAT_L001_R1.fastq P1-A01_GACGAT_L001_R2.fastq 10 output.fasta log.txt")
        sys.exit(1)
    
def checkArgumentMarker():
    # print(DnaUtils.reverseSequence(DnaUtils.complementDnaSequence("GATATTTTAGATATTTTCCGAAGGTACCATTTAA")))
    # sys.exit(1)
    if len(sys.argv) -1 != 4:
        print(sys.argv[0] + " [IUPACFile] [OligosFile] [output directory] [fasta file]")
        print("Exemple " + sys.argv[0] + " P1-A01_GACGAT_L001_R1.fastq P1-A01_GACGAT_L001_R2.fastq 10 output.fasta log.txt")
        sys.exit(1)
        
def checkArgumentMarkerReducer():
    if len(sys.argv) -1 != 6:
        print(sys.argv[0] + " [IUPACFile] [forward OligosFile] [reverse OligosFile] [marker directory] [marker name]")
        print("Exemple : python3 src/test.py data/iupac.csv data/amorce/oligos_f.oligos data/amorce/oligos_r.oligos data/marker/olon_COI olon_COI")
        for i in range(1, len(sys.argv)):
            print("[" + str(i) + "] " + sys.argv[i])
        sys.exit(1)
    
def mergeMap(globalmap, smallMap):
    for k in smallMap:
        if not k in globalmap:
            globalmap[k] = 0
        globalmap[k] += smallMap[k]

def generatePlot():
    checkArgumentMarkerReducer()
    logger.debug("IUPAC File : " +  sys.argv[2])
    iupacMap = DnaUtils.readIUPACFile(sys.argv[2])
    logger.debug("forward MarkerMap File : " +  sys.argv[3])
    logger.debug("reverse MarkerMap File : " +  sys.argv[4])
    forwardMarkerMap = DnaUtils.readOligosFile(sys.argv[3], iupacMap)
    reverseMarkerMap = DnaUtils.readOligosFile(sys.argv[4], iupacMap, True)
   
    rAllOutputFileName = sys.argv[6] + ".r.csv"
    fAllOutputFileName = sys.argv[6] + ".f.csv"
    
    fAll = {}
    rAll = {}
    
    for root, dirs, files in os.walk(sys.argv[5]):
        for f1le in files:
            filePath = os.path.join(sys.argv[5], f1le)
            logger.debug("file : " + filePath)
            f,r = DnaUtils.readmarker(filePath, forwardMarkerMap, reverseMarkerMap, sys.argv[6], False)
            
            mergeMap(rAll, r)
            mergeMap(fAll, f)
            
            indexOsPath = filePath.rfind(os.sep)
            if indexOsPath == -1:
                indexOsPath = 0
            else:
                indexOsPath += 1
            indexOfDot = filePath.rfind(".")
            if indexOfDot == -1:
                indexOfDot = len(filePath)
            rOutputFileName = filePath[indexOsPath:indexOfDot] + ".r.csv"
            fOutputFileName = filePath[indexOsPath:indexOfDot] + ".f.csv"
            DnaUtils.mapToCsv(r, rOutputFileName)
            DnaUtils.mapToCsv(f, fOutputFileName)
        
    DnaUtils.mapToCsv(rAll, rAllOutputFileName)
    DnaUtils.mapToCsv(fAll, fAllOutputFileName)

    
def testApplyThreshold():
    if len(sys.argv) -1 != 2:
        print(sys.argv[0] + " [fastq fileName] [threshold]")
        print("Exemple : python3 " + sys.argv[0] + " data/P1-A01_GACGAT_L001_R1.fastq 10")
        for i in range(1, len(sys.argv)):
            print("[" + str(i) + "] " + sys.argv[i])
        sys.exit(1)
    testQualityThreshold(sys.argv[1], int(sys.argv[2]))


def checkArgumentMarkerReducer():
    if len(sys.argv) -1 != 4:
        print(sys.argv[0] + " [IUPACFile] [OligosFile] [output directory] [fastq file]")
        for i in range(1, len(sys.argv)):
            print("[" + str(i) + "] " + sys.argv[i])
        sys.exit(1)
    

def readSingleFastQFile(fileName, threshold, logFileName, mapOligos, outputDir):
    f = open(fileName, "r")
    log = open(logFileName, "w")
    
    lineNumber = 0
    lastHeaderLineNumber    = 0
    lastSeparatorLineNumber = 0
    sequences  = []
    fastqSequence = None
    
    for line in f:
        line = line[:-1]
        if line[0] == "@" and lastSeparatorLineNumber != lineNumber - 1:
            if fastqSequence != None:
                #tryToFusion(fastqSequence, reverseFileName, threshold, log, seqOutput)
                fastqSequence.applyDrasticThreshold(threshold)
                DnaUtils.splitMarker1(mapOligos, outputDir, fastqSequence, fileName)
                
            fastqSequence = FastqSequence.readFasqSequenceHeader(line)
            lastHeaderLineNumber = lineNumber
            
            logging.debug(line)
        elif  lastHeaderLineNumber == lineNumber - 1:
            fastqSequence.setSequence(line)
        elif line[0] == "+":
            lastSeparatorLineNumber = lineNumber
        elif lastSeparatorLineNumber == lineNumber - 1:
            fastqSequence.setQuality(line)
            
        lineNumber += 1
    #tryToFusion(fastqSequence, reverseFileName, threshold, log, seqOutput)
    fastqSequence.applyDrasticThreshold(threshold)
    DnaUtils.splitMarker1(mapOligos, outputDir, fastqSequence, fileName)
    f.close()
    log.close()
    
    
def newPipeline():
    checkArgumentMarkerReducer()
    m = DnaUtils.readIUPACFile(sys.argv[1])
    mapOligos = DnaUtils.readOligosFile(sys.argv[2], m)
    readSingleFastQFile(sys.argv[4], 10, "log.txt", mapOligos, sys.argv[3])
    
if __name__ == '__main__':
    fileConfig('src/logging_config.ini')
    logger = logging.getLogger()
    newPipeline()
    #logging.basicConfig(filename='example.log',format='%(asctime)s:%(levelname)s:%(message)s', level=logging.CRITICAL)
    #checkArgument()
    #readFastQFile(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5])
    #print(DnaUtils.readIUPACFile(sys.argv[1]))
    # s fait du sort que le code sor
    #iupacMap = DnaUtils.readIUPACFile(sys.argv[1])
    #l = []
    #print(DnaUtils.rec_generate(sys.argv[2], 0, iupacMap, l));
    #sys.exit(1)
    
    # 21/03/2016
    #checkArgumentMarker()
    #m = DnaUtils.readIUPACFile(sys.argv[1])
    #mapOligos = DnaUtils.readOligosFile(sys.argv[2], m)
    #DnaUtils.splitMarker(mapOligos, sys.argv[3], sys.argv[4])
    
    #04/04/2016
    #generatePlot()
    #fastq = FastqSequence.readFasqSequenceHeader("@MISEQ3:7:000000000-MERGE3:1:1101:14204:1480 1:N:0:GTCTAC")
    #fastq.test()
    
