# coding: utf8
from DanaError import DanaError 
import DnaUtils

class FastqSequence:
    READ_DIRECTION_LEFT_TO_RIGHT = 1
    READ_DIRECTION_RIGHT_TO_LEFT = 2
    QUALITY_SCORE = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
    
    def __init__(self, idInstrumentName, runId, flowcellId, flowcellLane, tileNumber, xCoordinate, yCoordinate, readDirection, readFiltered, controlBit, sequenceIndex, sequence, quality):
        self.idInstrumentName = idInstrumentName
        self.runId          = runId
        self.flowcellId     = flowcellId
        self.flowcellLane   = flowcellLane
        self.tileNumber     = tileNumber
        self.xCoordinate    = xCoordinate
        self.yCoordinate    = yCoordinate
        self.readFiltered   = readFiltered
        self.setControlBit(controlBit)
        self.setReadDirection(readDirection)
        self.sequenceIndex  = sequenceIndex
        
        self.sequence       = sequence
        self.setQuality(quality)
        
    def setControlBit(self, value):
        if not isinstance(value, int):
            self.controlBit = int(value)
        
    def setReadDirection(self, value):
        readDirection = FastqSequence.READ_DIRECTION_LEFT_TO_RIGHT
        if isinstance(value, str):
            readDirection = int(value)
            
        if self.controlBit == 0:
            if readDirection != FastqSequence.READ_DIRECTION_LEFT_TO_RIGHT and readDirection != FastqSequence.READ_DIRECTION_RIGHT_TO_LEFT:
                raise DanaError("ReadDirection with invalid value : " + str(readDirection))
            else:
                self.readDirection = readDirection
        else:
            self.readDirection = readDirection
        
    def setQuality(self, value):
        if value == None:
            return
        if isinstance(value, str):
            if self.readDirection == FastqSequence.READ_DIRECTION_LEFT_TO_RIGHT:
                self.quality = qualityScoretoIntArray(value)
            else:
                reverse = value[::-1]
                self.quality = qualityScoretoIntArray(reverse)
        else:
            self.quality = value
        
    def setSequence(self, value):
        if self.readDirection == FastqSequence.READ_DIRECTION_LEFT_TO_RIGHT:
            self.sequence = value
        else:
            self.sequence = DnaUtils.complementAndReverseDnaSequence(value)
            
    def applyDrasticThreshold(self, threshold):
        try:
            index = self.quality.index(threshold)
            self.sequence = self.sequence[:index]
            del self.quality[index:]
        except ValueError:
            return
        
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
        string += str(self.readDirection)
        string += "\n"
        
        string += "readFiltered"
        string += " = "
        string += self.readFiltered
        string += "\n"
        
        string += "controlBit"
        string += " = "
        string += str(self.controlBit)
        string += "\n"
        
        string += "sequenceIndex"
        string += " = "
        string += self.sequenceIndex
        string += "\n"
        
        string += "sequence :\n"
        string += str(self.sequence)
        string += "\n"
        
        string += "quality :\n"
        string += str(self.quality)
        string += "\n"
        return string

def qualityScoretoIntArray(value):
    res = []
    for l in value:
        res.append(FastqSequence.QUALITY_SCORE.find(l))
    return res
    
def readFasqSequenceHeader(line):
    #print(line)
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
        
