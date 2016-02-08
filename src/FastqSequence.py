# coding: utf8

class FastqSequence:
    def __init__(self, idInstumentName, runId, flowcellId, tileNumber, xCoordinate, yCoordinate, readDirection, readFiltered, controlBit, sequenceIndex, sequence, quality):
        self.idInstrumentName = idInstrumentName
        self.runId          = runId
        self.flowcellId     = flowcellId
        self.tileNumber     = tileNumber
        self.xCoordinate    = xCoordinate
        self.yCoordinate    = yCoordinate
        self.readDirection  = readDirection
        self.readFiltered   = readFiltered
        self.controlBit     = controlBit
        self.sequenceIndex  = sequenceIndex
        self.sequence       = sequence
        self.quality        = quality   
    
def readFasqSequenceHeader(line):
    s = line.split(":")
    for e in s:
        print(e)
        