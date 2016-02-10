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
        else:
            self.controlBit = value
        
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
            if self.readDirection == FastqSequence.READ_DIRECTION_LEFT_TO_RIGHT:
                # index = self.quality.rindex(threshold)
                index = len(self.sequence) - 1
                while index >=0 and self.quality[index] <= threshold:
                    index -= 1
                index += 1
                # print("applyDrasticThreshold index : " + str(index))
                self.sequence = self.sequence[:index]
                del self.quality[index:]
            elif self.readDirection == FastqSequence.READ_DIRECTION_RIGHT_TO_LEFT:
                index = 0
                while index < len(self.quality) and self.quality[index] <= threshold:
                    index += 1
                # print(self.sequence)
                # print("applyDrasticThreshold reverse index : " + str(index))
                self.sequence = self.sequence[index:]
                # print(self.sequence)
                del self.quality[:index]
                
        except ValueError:
            return
    
    def getLineHeader(self):
        return "@" + self.idInstrumentName + ":" + self.runId + ":" + self.flowcellId + ":" + self.flowcellLane + ":" + self.tileNumber + ":" + self.xCoordinate + ":" + self.yCoordinate + " " + str(self.readDirection) + ":" + self.readFiltered + ":" + str(self.controlBit) + ":" + self.sequenceIndex
    
    def getReverseLineHeader(self):
        reverseReadDirection = FastqSequence.READ_DIRECTION_LEFT_TO_RIGHT
        if self.readDirection == FastqSequence.READ_DIRECTION_LEFT_TO_RIGHT:
            reverseReadDirection = FastqSequence.READ_DIRECTION_RIGHT_TO_LEFT
        return "@" + self.idInstrumentName + ":" + self.runId + ":" + self.flowcellId + ":" + self.flowcellLane + ":" + self.tileNumber + ":" + self.xCoordinate + ":" + self.yCoordinate + " " + str(reverseReadDirection) + ":" + self.readFiltered + ":" + str(self.controlBit) + ":" + self.sequenceIndex
      
    def __eq__(self, other):
        if isinstance(other, FastqSequence):
            return other.tileNumber == self.tileNumber and other.xCoordinate == self.xCoordinate and other.yCoordinate == self.yCoordinate
        return False
        
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

def fusion(fastq1, reverseFastq):
    fastqFusion = FastqSequence(fastq1.idInstrumentName, fastq1.runId, fastq1.flowcellId, fastq1.flowcellLane, fastq1.tileNumber, fastq1.xCoordinate, fastq1.yCoordinate, fastq1.readDirection, fastq1.readFiltered, fastq1.controlBit, fastq1.sequenceIndex, None, None)
    lastReverse10 = reverseFastq.sequence[0:10]
    found = fastq1.sequence.rfind(lastReverse10)
    sequence = fastq1.sequence[:found] + reverseFastq.sequence
    
    fastqFusion.setSequence(sequence)
    
    qualityList = []
    qualityList.extend(fastq1.quality[:found])
    qualityList.extend(reverseFastq.quality)
    fastqFusion.setQuality(qualityList)
    return fastqFusion

def matchFastqSequence(fastq1, fastq2):
    # print(fastq1.sequence)
    # print(fastq2.sequence)
    classicFastq = fastq1
    reverseFastq = fastq2
    if fastq2.readDirection == FastqSequence.READ_DIRECTION_RIGHT_TO_LEFT:
        reverseFastq = fastq2
        classicFastq = fastq1
    else:
        reverseFastq = fastq1
        classicFastq = fastq2
    
    numberChar = 10
    lastReverse10 = reverseFastq.sequence[0:numberChar]
    found = classicFastq.sequence.rfind(lastReverse10)
    if found != -1:
        resteAVerifier = len(classicFastq.sequence) - found - numberChar
        # print("resteAVerifier " + str(resteAVerifier) + " caracteres")
        # print(" ->" + classicFastq.sequence[found + 10:])
        j = numberChar
        match = True
        for i in range(found + numberChar, len(classicFastq.sequence)):
            # print("[" + str(i) + "] : " + classicFastq.sequence[i] + "/" + reverseFastq.sequence[j])
            if classicFastq.sequence[i] != reverseFastq.sequence[j]:
                match = False
                break
            j += 1
        if match:
            return fusion(classicFastq, reverseFastq)
    return None

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
        
