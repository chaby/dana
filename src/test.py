# coding: utf8
import DnaUtils
import FastqSequence
import logging
import sys
import re

def tryToFusion(fastqSequence, reverseFileName, threshold, log, seqOutput):
    fastqSequence
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
    
if __name__ == '__main__':
    logging.basicConfig(filename='example.log',format='%(asctime)s:%(levelname)s:%(message)s', level=logging.CRITICAL)
    checkArgument()
    readFastQFile(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5])
    #fastqSequence = FastqSequence.readFasqSequenceHeader("@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG")
    #print(fastqSequence)
