# coding: utf8
import FastqSequence
import sys

def readFastQFile(fileName):
    f = open(fileName, "r")
    
    for line in f:
        line = line[:-1]
        if line[0] == "@":
            FastqSequence.readFasqSequenceHeader(line)
    f.close()
    
if __name__ == '__main__':
    readFastQFile(sys.argv[1])