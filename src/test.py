# coding: utf8
import FastqSequence
import sys
import re

def tryToFusion(fastqSequence, reverseFileName, threshold):
    fastqSequence
    f = open(reverseFileName, "r")
    reverseLineHeader = fastqSequence.getReverseLineHeader()
    
    lastHeaderLineNumber    = 0
    lastSeparatorLineNumber = 0
    parseFastq     = False
    rFastqSequence = None
    lineNumber     = 0
    
    for line in f:
        line = line[:-1]
        lineNumber += 1
        
        if line[0] == "@" and lastSeparatorLineNumber != lineNumber - 1:
            if parseFastq:
                rFastqSequence.applyDrasticThreshold(threshold)
                s = FastqSequence.matchFastqSequence(fastqSequence, rFastqSequence)
                if s != None:
                    print(s)
                else:
                    print("Can't merge " + fastqSequence.getLineHeader() + " with " + rFastqSequence.getLineHeader())
                parseFastq = False
                
            if line == reverseLineHeader:
                parseFastq = True
                
            if parseFastq:
                rFastqSequence = FastqSequence.readFasqSequenceHeader(line)
                lastHeaderLineNumber = lineNumber
        elif parseFastq:
            if  lastHeaderLineNumber == lineNumber - 1:
                rFastqSequence.setSequence(line)
            elif line[0] == "+":
                lastSeparatorLineNumber = lineNumber
            elif lastSeparatorLineNumber == lineNumber - 1:
                rFastqSequence.setQuality(line)                
            
    f.close()

def readFastQFile(fileName, reverseFileName, threshold):
    f = open(fileName, "r")
    
    lineNumber = 0
    lastHeaderLineNumber    = 0
    lastSeparatorLineNumber = 0
    sequences  = []
    fastqSequence = None
    
    for line in f:
        line = line[:-1]
        if line[0] == "@" and lastSeparatorLineNumber != lineNumber - 1:
            if fastqSequence != None:
                tryToFusion(fastqSequence, reverseFileName, threshold)
                
            fastqSequence = FastqSequence.readFasqSequenceHeader(line)
            lastHeaderLineNumber = lineNumber
        elif  lastHeaderLineNumber == lineNumber - 1:
            fastqSequence.setSequence(line)
        elif line[0] == "+":
            lastSeparatorLineNumber = lineNumber
        elif lastSeparatorLineNumber == lineNumber - 1:
            fastqSequence.setQuality(line)
            
        lineNumber += 1
    tryToFusion(fastqSequence, reverseFileName, threshold)
    f.close()
    
    #for s in sequences:
    #    s.applyDrasticThreshold(threshold)
    
    #f = FastqSequence.matchFastqSequence(sequences[0], sequences[1])
    #print(f.sequence)
    #print(sequences[0] == sequences[1])

#def searchSequence(fastqSequence, fileName):
    
def checkArgument():
    
    if len(sys.argv) -1 != 3:
        print(sys.argv[0] + " [R1.fastq] [R2.fastq] [drastic threshold]")
        print("Exemple " + sys.argv[0] + " P1-A01_GACGAT_L001_R1.fastq P1-A01_GACGAT_L001_R2.fastq 10")
        sys.exit(1)
    
if __name__ == '__main__':
    checkArgument()
    readFastQFile(sys.argv[1], sys.argv[2], int(sys.argv[3]))
    #fastqSequence = FastqSequence.readFasqSequenceHeader("@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG")
    #print(fastqSequence)
