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
                logger.debug("before threshold " + fastqSequence.sequence)
                fastqSequence.applyDrasticThreshold(threshold)
                logger.debug("after threshold " + fastqSequence.sequence)
                logger.debug("after threshold " + DnaUtils.complementAndReverseDnaSequence(fastqSequence.sequence))
                DnaUtils.splitMarker1(mapOligos, outputDir, fastqSequence, fileName)
                
            fastqSequence = FastqSequence.readFasqSequenceHeader(line)
            lastHeaderLineNumber = lineNumber
            
            logging.debug(line)
        elif  lastHeaderLineNumber == lineNumber - 1:
            logging.debug("seq in file : " + line)
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
    # for k in mapOligos:
    #     for e in mapOligos[k]:
    #         print(e)
    # sys.exit(1)
    readSingleFastQFile(sys.argv[4], 10, "log.txt", mapOligos, sys.argv[3])
    
def readMultipleThesholdWell(fileName):
    f = open(fileName, "r")
    m = {}
    g = {}
    forward  = 0
    reverse  = 1
    lineNumber = 0
    currentWellName = None
    for line in f:
        lineNumber += 1
        line = line[:-1]
        print(line)
        a = line.split("\t")
        markerName = a[0]
        direction  = a[1]
        theshold   = int(a[2])
        if not markerName in m:
            m[markerName] = []
            m[markerName].append(-1)
            m[markerName].append(-1)
            
        if direction == "forward":
            m[markerName][forward] = theshold
        elif direction == "reverse":
            m[markerName][reverse] = theshold
            newName = markerName + "_" + str(m[markerName][forward]) + "_" + str(m[markerName][reverse])
            if not markerName in g:
                g[markerName] = {}
            g[markerName][newName] = m[markerName][:]
        else:
            logger.error("direction : \"" + "\" invalid valud in file : " + fileName + " at line [" + str(lineNumber) + "]")
    f.close()
    return g
    
def readThesholdWell(fileName):
    f = open(fileName, "r")
    m = {}
    forward  = 0
    reverse  = 1
    lineNumber = 0
    for line in f:
        lineNumber += 1
        line = line[:-1]
        print(line)
        a = line.split("\t")
        markerName = a[0]
        direction  = a[1]
        theshold   = int(a[2])
        if not markerName in m:
            m[markerName] = []
            m[markerName].append(-1)
            m[markerName].append(-1)
            
        if direction == "forward":
            m[markerName][forward] = theshold
        elif direction == "reverse":
            m[markerName][reverse] = theshold
        else:
            logger.error("direction : \"" + "\" invalid valud in file : " + fileName + " at line [" + str(lineNumber) + "]")
    f.close()
    return m

def checkArgumentEcorAndOlon():
    if len(sys.argv) -1 != 4:
        print(sys.argv[0] + " [well forward] [well reverse] [threshold well] [well output]")
        for i in range(1, len(sys.argv)):
            print("[" + str(i) + "] " + sys.argv[i])
        sys.exit(1)
        
def searchCorrespondingReverseSeq(header, reverseFilehandler, wellThresholds):
    reverseFilehandler.seek(0)
    
    reverseHeader = header.replace(" 1:", " 2:")
    lastColon     = reverseHeader.rfind(":")
    reverseHeader = reverseHeader[:lastColon]
    
    #logger.debug("[searchCorrespondingReverseSeq] reverseHeader : " + reverseHeader)
    header = None
    data   = ""
    for line in reverseFilehandler:
        line = line[:-1]
        if len(line) > 0 and line[0] == "@":
            if line.find(reverseHeader) == 0:
                header = line
            # last header was the good one
            # i hava fish to read fasta
            elif header != None:
                break
            else:
                data = ""
        else:
            data += line
    
    return (header, data)
       
            
def cutFusionAndMerge(headerF, dataF, headerR, dataR, wellThresholds, log, outputFile):
    logger.debug("[cutFusionAndMerge] " + headerF + ", " + headerR + ", " + str(wellThresholds))
    
    forward  = 0
    reverse  = 1
    
    dataF = dataF[0:wellThresholds[forward]]
    #a = "ATAACGCTGTTATCCCTGCGGTAACTTGTTCTTTTGATCACTGTAAGTGGATCACACCTTCATTTTTATGATTTAAGAAAAACAATTCTTTTATTTTAGGTTAATATAACCATATAGTAGCGGAGGATTTTCTTTCTCCGGGATTGCCCCAATCAAAGCTTGTTTCAATTTGCCATGCTCTAGGCCTACTATTTCTATTATATTAGTTAGGGCTAATAGTAAATAACAATTAAAATTCAACTACAGCTCG"
    #dataR = DnaUtils.complementAndReverseDnaSequence(dataR[0:wellThresholds[reverse]])
    cdataR = DnaUtils.complementAndReverseDnaSequence(dataR)
    cdataR = cdataR[:wellThresholds[reverse]]
    #cdataR = DnaUtils.complementAndReverseDnaSequence(cdataR)
    #data   = merge(dataF, cdataR)
    # print(dataF)
    #print(DnaUtils.complementAndReverseDnaSequence(cdataR))
    
    fastqSequenceForward = FastqSequence.createFasqSequenceHeader(headerF)
    fastqSequenceForward.setSequence(dataF, True)
    fastqSequenceForward.quality = None
    
    fastqSequenceReverse = FastqSequence.createFasqSequenceHeader(headerR)
    fastqSequenceReverse.setSequence(cdataR)
    fastqSequenceReverse.quality = None
    
    newFastQ = FastqSequence.matchFastqSequence(fastqSequenceForward, fastqSequenceReverse)
    
    logging.debug(newFastQ.getLineHeader() + "\t" + str(newFastQ.typeFusion))
    if newFastQ.typeFusion == newFastQ.TYPE_FUSION_OK:
        log.write("Merge " + fastqSequenceForward.getLineHeader() + " with " + fastqSequenceReverse.getLineHeader() + "\n")
    elif newFastQ.typeFusion == newFastQ.TYPE_FUSION_NO_MATCH:
        log.write("Forced merge " + fastqSequenceForward.getLineHeader() + " with " + fastqSequenceReverse.getLineHeader() + "\n")
    elif newFastQ.typeFusion == newFastQ.TYPE_FUSION_PARTIAL:
        log.write("Partial merge " + fastqSequenceForward.getLineHeader() + " with " + fastqSequenceReverse.getLineHeader() + "\n")
    else:
        log.write("Can't merge " + fastqSequenceForward.getLineHeader() + " with " + fastqSequenceReverse.getLineHeader() + "\n")
        
    #outputFile.write(">" + newFastQ.getLineHeader() + "\n")
    #outputFile.write(newFastQ.sequence + "\n")
    #sys.exit(1)
    return newFastQ
    
    
def tryToFusionMarkerFile(forwardMarkerFile, reverseMarkerFile, wellThresholds, log, logDoublon, outputFile, fusionMap):
    logger.debug("[tryToFusionMarkerFile] " + forwardMarkerFile + ", " + reverseMarkerFile + ", " + str(wellThresholds))
    f = open(forwardMarkerFile, "r")
    r = open(reverseMarkerFile, "r")
    
    headerF = None
    dataF   = ""
    headerR = None
    dataR   = None
        
    for line in f:
        line = line[:-1]
        if len(line) > 0 and line[0] == "@":
            if headerF == None:
                headerF = line
            elif line != headerF:
                (headerR, dataR) = searchCorrespondingReverseSeq(headerF, r, wellThresholds)
                #if headerR == None:
                #    logger.error("[searchCorrespondingReverseSeq] " + headerF + " not found in " + reverseMarkerFile)
                #else:
                if headerR != None:
                    fastq = cutFusionAndMerge(headerF, dataF, headerR, dataR, wellThresholds, log, outputFile)
                    appendToDoublonMap(fastq, fusionMap)
                        
                headerF = line
                dataF   = ""
        else:
            dataF += line
            
    (headerR, dataR) = searchCorrespondingReverseSeq(headerF, r, wellThresholds)
    if headerR == None:
        logger.error("[searchCorrespondingReverseSeq] END " + headerF + " not found in " + reverseMarkerFile)
    else:
        fastq = cutFusionAndMerge(headerF, dataF, headerR, dataR, wellThresholds, log, outputFile)
    appendToDoublonMap(fastq, fusionMap)
    
    r.close()
    f.close()
    
def appendToDoublonMap(fastq, fusionMap):
    if fastq.sequence in fusionMap:
        fusionMap[fastq.sequence][1] = fusionMap[fastq.sequence][1] +1
        fusionMap[fastq.sequence][2].append(fastq.getLineHeader())
    else:
        fusionMap[fastq.sequence] = []
        fusionMap[fastq.sequence].append(fastq)
        fusionMap[fastq.sequence].append(1)
        emptyDoublon = []
        fusionMap[fastq.sequence].append(emptyDoublon)
    
        
def tryToFusionWell(wellForward, wellReverse, wellThresholds, outputDir):
    logger.info("[tryToFusionWell] " + wellForward + ", " + wellReverse + ", " + str(wellThresholds))
    for root, dirs, files in os.walk(wellForward):
        for f in files:
            logger.debug("[tryToFusionWell] outputFile : " + os.path.join(outputDir, f))
            logMergeFileName   = os.path.join(outputDir, f[:f.rfind(".")] + ".merge.log")
            logDoublonFileName = os.path.join(outputDir, f[:f.rfind(".")] + ".doublon.log")
            logCountFileName = os.path.join(outputDir, f[:f.rfind(".")] + ".count.log")
            m = {}
            logger.debug("[tryToFusionWell] logFile    : " + logMergeFileName)
            of = open(os.path.join(outputDir, f), "w")
            
            logMerge = open(logMergeFileName, "w")
            logDoublon = open(logDoublonFileName, "w")
            logCount = open(logCountFileName, "w")
            # logger.debug("[tryToFusionWell] file : " + f + " " + os.path.join(wellReverse, f) + " : " + str(os.path.exists(os.path.join(wellReverse, f))))
            if (os.path.exists(os.path.join(wellReverse, f))):
                tryToFusionMarkerFile(os.path.join(wellForward, f), os.path.join(wellReverse, f), wellThresholds, logMerge, logDoublon, of, m)
            
            for k in m:
                a = m[k]
                fastq = a[0]
                count = a[1]
                elms  = a[2]
                if count == 1:
                    logDoublon.write(fastq.getLineHeader() + "\t1\tdelete\n");
                elif count > 1:
                    of.write(">" + fastq.getLineHeader()+"\n")
                    of.write(fastq.sequence + "\n")
                    logCount.write(fastq.getLineHeader() + "\t"+str(count)+"\tkeep\n");
                    for elm in elms:
                        logDoublon.write(elm + "\t" + str(count) + "\tidentical\t" + fastq.getLineHeader() + "\n")
            
            logCount.close()
            logDoublon.close()
            logMerge.close()
            of.close()
    
    
def testEcorAndOlon():
    checkArgumentEcorAndOlon()
    wellForwardDirectory = sys.argv[1]
    wellReverseDirectory = sys.argv[2]
    #wellThresholdMap     = readThesholdWell(sys.argv[3])
    wellThresholdMap     = readMultipleThesholdWell(sys.argv[3])
    weelOutput           = sys.argv[4]
    
    #print(wellThresholdMap)
    for root, dirs, files in os.walk(wellForwardDirectory):
        for directory in dirs:
            if directory in wellThresholdMap:
                logger.debug("diractory in " + wellForwardDirectory + " : " + directory)
                if not os.path.exists(os.path.join(weelOutput, directory)):
                    os.makedirs(os.path.join(weelOutput, directory))
                    
                if (os.path.exists(os.path.join(wellReverseDirectory, directory))):
                    logger.debug("reverse diractory in : " + os.path.join(wellReverseDirectory, directory))
                    if len(wellThresholdMap[directory]) == 1:
                        wellName, wellValues = wellThresholdMap[directory].popitem()
                        
                        tryToFusionWell(os.path.join(wellForwardDirectory, directory), os.path.join(wellReverseDirectory, directory), wellValues, os.path.join(weelOutput, directory))
                    else:
                        for wells in wellThresholdMap[directory]:
                            if not os.path.exists(os.path.join(weelOutput, directory, wells)):
                                os.makedirs(os.path.join(weelOutput, directory, wells))
                            tryToFusionWell(os.path.join(wellForwardDirectory, directory), os.path.join(wellReverseDirectory, directory), wellThresholdMap[directory][wells], os.path.join(weelOutput, directory, wells))
    
if __name__ == '__main__':
    fileConfig('src/logging_config.ini')
    logger = logging.getLogger()
    #newPipeline()
    testEcorAndOlon()
    
    #print(DnaUtils.complementAndReverseDnaSequence("CGAGCTGTATCTGAATTTTAATTGTTATTTACTATTAGCCCTAACTAATATGATAGAAATAGTAGGCCTAAAGCATGGCAAATTGAAACAAGCTTTGATTGGGGCAATCCCGGAGAAAGAAAATCCTCCGCTACTATATGGTTATATTAACCTAAAATAAAAGAATTGTTTTTCTTAAATCATAAAAATGAAGGTGTGATCCACTTACAGTGATCAAAAGAACAAGTT"))
    
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
    
