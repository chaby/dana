# coding: utf8
from DanaError import DanaError 
import DnaUtils
import logging
import sys

class FastqSequence(object):
    READ_DIRECTION_LEFT_TO_RIGHT = 1
    READ_DIRECTION_RIGHT_TO_LEFT = 2
    
    TYPE_FUSION_OK = 0
    TYPE_FUSION_NO_MATCH = 1
    TYPE_FUSION_PARTIAL = 2
    
    QUALITY_SCORE = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~"
    
    def __init__(self, idInstrumentName, runId, flowcellId, flowcellLane, tileNumber, xCoordinate, yCoordinate, readDirection, readFiltered, controlBit, sequenceIndex, sequence, quality):
        self.logger = logging.getLogger(__name__)
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
        self.typeFusion     = -1
        
    def test(self):
        self.logger.debug("test debug")
        self.logger.info("test debug")
    
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
        
    def setSequence(self, value, forced=False):
        if forced:
            self.sequence = value
        else:
            if self.readDirection == FastqSequence.READ_DIRECTION_LEFT_TO_RIGHT:
                self.sequence = value
            else:
                self.sequence = DnaUtils.complementAndReverseDnaSequence(value)
        self.logger.debug("self.sequence : " + self.sequence)
                
    def isReverse(self):
        return self.readDirection == FastqSequence.READ_DIRECTION_RIGHT_TO_LEFT
        
    def applyDrasticThreshold(self, threshold):
        try:
            if self.readDirection == FastqSequence.READ_DIRECTION_LEFT_TO_RIGHT:
                
                index = len(self.sequence) - 1
                while index >=0 and self.quality[index] <= threshold:
                    index -= 1
                index += 1
                
                self.sequence = self.sequence[:index]
                del self.quality[index:]
            elif self.readDirection == FastqSequence.READ_DIRECTION_RIGHT_TO_LEFT:
                index = 0
                while index < len(self.quality) and self.quality[index] <= threshold:
                    index += 1
                                
                self.sequence = self.sequence[index:]
                self.logger.debug("[applyDrasticThreshold] REVERSE delete from 0 to " + str(index))
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

def partialFusion(fastq1, reverseFastq, size):
    fastqFusion = FastqSequence(fastq1.idInstrumentName, fastq1.runId, fastq1.flowcellId, fastq1.flowcellLane, fastq1.tileNumber, fastq1.xCoordinate, fastq1.yCoordinate, fastq1.readDirection, fastq1.readFiltered, fastq1.controlBit, fastq1.sequenceIndex, None, None)
    #sequence = fastq1.sequence + "XXXX" + reverseFastq.sequence
    sequence = fastq1.sequence
    for i in range(0,size):
        sequence += "N"
    sequence += reverseFastq.sequence[size:]
    fastqFusion.setSequence(sequence)
    qualityXXXX  = []
    for i in range(0,size):
        qualityXXXX.append(0)
        
    if fastq1.quality != None:
        qualityList = []
        qualityList.extend(fastq1.quality)
        qualityList.extend(qualityXXXX)
        qualityList.extend(reverseFastq.quality[size:])
        fastqFusion.setQuality(qualityList)
    else:
        fastqFusion.quality = None
    fastqFusion.typeFusion = FastqSequence.TYPE_FUSION_PARTIAL
    return fastqFusion
    
def noMatchFusion(fastq1, reverseFastq):
    fastqFusion = FastqSequence(fastq1.idInstrumentName, fastq1.runId, fastq1.flowcellId, fastq1.flowcellLane, fastq1.tileNumber, fastq1.xCoordinate, fastq1.yCoordinate, fastq1.readDirection, fastq1.readFiltered, fastq1.controlBit, fastq1.sequenceIndex, None, None)
    sequence = fastq1.sequence + "XXXX" + reverseFastq.sequence
    fastqFusion.setSequence(sequence)
    if fastq1.quality != None:
        qualityXXXX  = []
        qualityXXXX.append(0)
        qualityXXXX.append(0)
        qualityXXXX.append(0)
        qualityXXXX.append(0)
        qualityList = []
        qualityList.extend(fastq1.quality)
        qualityList.extend(qualityXXXX)
        qualityList.extend(reverseFastq.quality)
        fastqFusion.setQuality(qualityList)
    else:
        fastqFusion.quality = None
    fastqFusion.typeFusion = FastqSequence.TYPE_FUSION_NO_MATCH
    return fastqFusion

def fusion(fastq1, reverseFastq):
    fastqFusion = FastqSequence(fastq1.idInstrumentName, fastq1.runId, fastq1.flowcellId, fastq1.flowcellLane, fastq1.tileNumber, fastq1.xCoordinate, fastq1.yCoordinate, fastq1.readDirection, fastq1.readFiltered, fastq1.controlBit, fastq1.sequenceIndex, None, None)
    lastReverse10 = reverseFastq.sequence[0:10]
    found = fastq1.sequence.rfind(lastReverse10)
    sequence = fastq1.sequence[:found] + reverseFastq.sequence
    
    fastqFusion.setSequence(sequence)
    
    if fastq1.quality != None:
        qualityList = []
        qualityList.extend(fastq1.quality[:found])
        qualityList.extend(reverseFastq.quality)
        fastqFusion.setQuality(qualityList)
    else:
        fastqFusion.quality = None
    fastqFusion.typeFusion = FastqSequence.TYPE_FUSION_OK
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
        
    logging.debug(classicFastq.sequence)
    logging.debug(reverseFastq.sequence)
    
    numberChar = 10
    lastReverse10 = reverseFastq.sequence[0:numberChar]
    found = classicFastq.sequence.rfind(lastReverse10)
    logging.debug("lastReverse10 : " + lastReverse10)
    logging.debug("found : " + str(found))
    if found == -1:
        numberChar -=1
        while (numberChar) > 1:
            lastReverse10 = reverseFastq.sequence[0:numberChar]
            found = classicFastq.sequence.rfind(lastReverse10)
            logging.debug("lastReverse" + str(numberChar) + " : " + lastReverse10)
            logging.debug("found : " + str(found))
            if found != -1:
                resteAVerifier = len(classicFastq.sequence) - found - numberChar
                logging.debug("resteAVerifier " + str(resteAVerifier) + " caracteres")
                logging.debug(" ->" + classicFastq.sequence[found + 10:])
                allMatch = True
                mF = ""
                mR = ""
                for k in range(0, resteAVerifier):
                    if resteAVerifier < len(reverseFastq.sequence):
                        mF += classicFastq.sequence[found + k]
                        mR += reverseFastq.sequence[k]
                        if reverseFastq.sequence[k] != classicFastq.sequence[found + k]:
                            allMatch = False
                            break
                
                #print(classicFastq.sequence)
                #k = 0
                #s = ""
                #for k in range(0, found):
                #    s += "#"
                #print(s + reverseFastq.sequence)
                
                #print(mF)
                #print(mR)
                if not allMatch:
                    return noMatchFusion(classicFastq, reverseFastq)
                
                        
                
                return partialFusion(classicFastq, reverseFastq, numberChar)
                break
                
            numberChar -=1
            
    if found != -1:
        resteAVerifier = len(classicFastq.sequence) - found - numberChar
        # print("resteAVerifier " + str(resteAVerifier) + " caracteres")
        # print(" ->" + classicFastq.sequence[found + 10:])
        logging.debug("resteAVerifier " + str(resteAVerifier) + " caracteres")
        logging.debug(" ->" + classicFastq.sequence[found + 10:])
        j = numberChar
        match = True
        for i in range(found + numberChar, len(classicFastq.sequence)):
            # print("[" + str(i)  + "/" + str(len(classicFastq.sequence))+ "] : " + "[" + str(j)  + "/" + str(len(reverseFastq.sequence))+ "] : " + classicFastq.sequence[i] + "/" + reverseFastq.sequence[j])
            if j >= len(reverseFastq.sequence) or i > len(classicFastq.sequence) or classicFastq.sequence[i] != reverseFastq.sequence[j]:
                match = False
                break
            j += 1
        if match:
            return fusion(classicFastq, reverseFastq)
    logging.debug("No match")
    return noMatchFusion(classicFastq, reverseFastq)

def qualityScoretoIntArray(value):
    res = []
    for l in value:
        res.append(FastqSequence.QUALITY_SCORE.find(l))
    return res
    
def createFasqSequenceHeader(line):
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
        
