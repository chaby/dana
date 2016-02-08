# coding: utf8

class FastqSequence:
    def __init__(self, idInstrumentName, runId, flowcellId, flowcellLane, tileNumber, xCoordinate, yCoordinate, readDirection, readFiltered, controlBit, sequenceIndex, sequence, quality):
        self.idInstrumentName = idInstrumentName
        self.runId          = runId
        self.flowcellId     = flowcellId
        self.flowcellLane   = flowcellLane
        self.tileNumber     = tileNumber
        self.xCoordinate    = xCoordinate
        self.yCoordinate    = yCoordinate
        self.readDirection  = readDirection
        self.readFiltered   = readFiltered
        self.controlBit     = controlBit
        self.sequenceIndex  = sequenceIndex
        
        self.sequence       = sequence
        self.quality        = quality
        
    def __str__(self):
        string = ""
        string += "idInstrumentName"
        string += " = "
        string += self.idInstrumentName
        string += "\n"
        
        string += "runId"
        string += " = "
        string += self.runId
        string += "\n"
        
        string += "flowcellId"
        string += " = "
        string += self.flowcellId
        string += "\n"
        
        string += "flowcellLane"
        string += " = "
        string += self.flowcellLane
        string += "\n"
        
        string += "tileNumber"
        string += " = "
        string += self.tileNumber
        string += "\n"
        
        string += "xCoordinate"
        string += " = "
        string += self.xCoordinate
        string += "\n"
        
        string += "yCoordinate"
        string += " = "
        string += self.yCoordinate
        string += "\n"
        
        string += "readDirection"
        string += " = "
        string += self.readDirection
        string += "\n"
        
        string += "readFiltered"
        string += " = "
        string += self.readFiltered
        string += "\n"
        
        string += "controlBit"
        string += " = "
        string += self.controlBit
        string += "\n"
        
        string += "sequenceIndex"
        string += " = "
        string += self.sequenceIndex
        string += "\n"
        return string
    
def readFasqSequenceHeader(line):
    print(line)
    s = line.split(":")
    idInstrumentName = s[0][1:]
    runId            = s[1]
    flowcellId       = s[2]
    flowcellLane     = s[3]
    tileNumber       = s[4]
    xCoordinate      = s[5]
    ss = s[6].split(" ")
    yCoordinate      = ss[0]
    readDirection    = ss[1]
    readFiltered     = s[7]
    controlBit       = s[8]
    sequenceIndex    = s[9]
    fastqSequence = FastqSequence(idInstrumentName, runId, flowcellId, flowcellLane, tileNumber, xCoordinate, yCoordinate, readDirection, readFiltered, controlBit, sequenceIndex, None, None)
    return fastqSequence
        
