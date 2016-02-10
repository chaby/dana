# coding: utf8
import FastqSequence
import sys

def readFastQFile(fileName, threshold):
    f = open(fileName, "r")
    
    lineNumber = 0
    lastHeaderLineNumber    = 0
    lastSeparatorLineNumber = 0
    sequences  = []
    fastqSequence = None
    
    for line in f:
        line = line[:-1]
        #print(str(lineNumber) + "\t" + str(lastHeaderLineNumber) + "\t" + line)
        if line[0] == "@":
            if fastqSequence != None:
                sequences.append(fastqSequence)
            fastqSequence = FastqSequence.readFasqSequenceHeader(line)
            lastHeaderLineNumber = lineNumber
        elif  lastHeaderLineNumber == lineNumber - 1:
            fastqSequence.setSequence(line)
        elif line[0] == "+":
            lastSeparatorLineNumber = lineNumber
        elif lastSeparatorLineNumber == lineNumber - 1:
            fastqSequence.setQuality(line)
            
        lineNumber += 1
    sequences.append(fastqSequence)
    f.close()
    
    for s in sequences:
        s.applyDrasticThreshold(threshold)
    
    f = FastqSequence.matchFastqSequence(sequences[0], sequences[1])
    print(f.sequence)
    print(sequences[0] == sequences[1])

#def searchSequence(fastqSequence, fileName):
    
    
    
if __name__ == '__main__':
    readFastQFile(sys.argv[1], int(sys.argv[2]))
    #fastqSequence = FastqSequence.readFasqSequenceHeader("@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG")
    #print(fastqSequence)
