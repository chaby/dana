# coding: utf8
import os
import os.path

def complementDnaSequence(sequence):
    complement = ""
    for l in sequence:
        if l == "A":
            complement += "T"
        elif l == "T":
            complement += "A"
        elif l == "C":
            complement += "G"
        elif l == "G":
            complement += "C"
        else:
            complement += l
    return complement
    

    
def complementAndReverseDnaSequence(sequence):
    return complementDnaSequence(sequence)[::-1]

def reverseSequence(sequence):
    return sequence[::-1]
    
def readIUPACFile(fileName):
    uf = open(fileName, "r")
    iupacMap = {}
    for line in uf:
        line = line[:-1]
        
        s = line.split("\t")
        correspondance = s[1].split(" or ")
        iupacMap[s[0]] = correspondance
        #print(line)
    uf.close()
    return iupacMap

#def generateIUPACPossibilities(string, iupacMap):
    #for c in string:
        #if c in iupacMap:
            
def rec_generate(string, pos, iupacMap, listz):
    #print("[rec_generate] " + string + "/" + str(pos) + "/" + str(len(string)) + " / " + str(pos > len(string)) + " -> " + str(listz))
    
    if pos >= len(string):
        return listz
        
    #if len(listz) == 0:
    #    listz.append("")
        
    if string[pos] in iupacMap:
        #print("string[pos] = \"" + string[pos] + " \" in iupacMap len(listz) : " + str(len(listz)))
        
        oldSize = len(listz)
        if oldSize == 0:
            for elm in iupacMap[string[pos]]:
                listz.append(elm)
        else:
            ### pour chaque element de la liste on rajoute un element pour chaque
            ### possibilite de IUPAC
            for k in range(0, len(iupacMap[string[pos]]) - 1 ):
                for j in range(0,oldSize):
                    listz.append(listz[j])
            #print("# string[pos] in iupacMap " + str(listz) + "/" + str(len(iupacMap[string[pos]])))
        
            ### each IUPAC possibilities, add it at position modulo the old list size
            for l in range(0, len(iupacMap[string[pos]])):
                for i in range(0, oldSize):
                    index = l * oldSize + i
                    
                    #print("i:" + str(i) +" => " + str(index) + " ajoute " + iupacMap[string[pos]][l])
                    #print("Ajoute a " + listz[i * oldSize + l] + " la lettre " + iupacMap[string[pos]][l])
                    listz[index] = listz[index] + iupacMap[string[pos]][l]
                
        #print(listz)
    else:
        #print("string[pos] NOT in iupacMap")
        if len(listz) == 0:
            listz.append("")
        i = 0
        for elm in listz:
            listz[i] = elm + string[pos]
            #print("elm : " + listz[i])
            i += 1
        #print(listz)
    return rec_generate(string, pos + 1, iupacMap, listz)
    
            
    
def readOligosFile(fileName, iupacMap, reverse=False):
    uf = open(fileName, "r")
    m = {}
    for line in uf:
        line = line[:-1]
        
        s = line.split("\t")
        marker = s[0].upper()
        if reverse:
            marker = complementAndReverseDnaSequence(marker)
            # print ("reverse : " + marker + "  " + str(rec_generate(marker, 0, iupacMap, [])))
        m[s[1]] = rec_generate(marker, 0, iupacMap, [])
        
    uf.close()
    return m
    
def getWellName(string):
    lastIndexOfSlash = string.rfind("/")
    if lastIndexOfSlash == -1:
        lastIndexOfSlash = 0
    firstIndexOfUnderscore = string.find("_", lastIndexOfSlash)
    return string[lastIndexOfSlash + 1:firstIndexOfUnderscore]
    
def splitMarker(markerMap, outputDirectory, fileName):
    wellFileName = getWellName(fileName)
    
    for directory in markerMap.keys():
        if not os.path.exists(os.path.join(outputDirectory,directory)):
            os.makedirs(os.path.join(outputDirectory,directory))
        fastaFile = open(fileName, "r")
        header = ""
        for line in fastaFile:
            line = line[:-1]
            if line.find(">") == 0:
                header = line
            else:
                for elm in markerMap[directory]:
                    index = line.find(elm)
                    if index != -1 and index < 50:
                        #print (elm + "("+str(index)+") trouve dans " + header + " ajout dans " + os.path.join(outputDirectory,directory,wellFileName) + ".fasta")
                        f = open(os.path.join(outputDirectory,directory,wellFileName) + ".fasta", "a")
                        f.write(header + os.linesep)
                        f.write(line + os.linesep)
                        f.close()
        fastaFile.close()

# TODO : que faire si on ne trouve pas la sequence marker et pour reverse
def readmarker(markerFile, forwardMarkerMap, reverseMarkerMap, markerName, trace=True):
    NMap = {}
    RNMap = {}
    NMap[-1] = 0
    RNMap[-1] = 0
    mf = open(markerFile, "r")
    for line in mf:
        line = line[:-1]
        if line.find(">") == -1:
            
            forwardMarkerList = forwardMarkerMap[markerName]
            reverseMarkerList = reverseMarkerMap[markerName]
            found = False
            for elm in forwardMarkerList:
                if line.find(elm) == 0:
                    found = True
                    break
            # delete the first characters corresponding to the marker sequence
            line = line[len(forwardMarkerList[0]):]
            
            
            found = False
            for elm in reverseMarkerList:
                index = line.find(elm)
                if index != -1 and len(elm) + index == len(line):
                    found = True
                    break
            
            # delete the last characters corresponding to the reverse marker sequence
            line = line[0:len(line) - len(reverseMarkerList[0])]
            if trace:
                print(line)
            
            indexOfN = line.find("N")
            indexOfX = line.find("X")
            
            if not indexOfN in NMap:
                NMap[indexOfN] = 0
                
            if not indexOfX in NMap:
                NMap[indexOfX] = 0
                
            if indexOfN == -1 and indexOfX == -1:
                NMap[-1] += 1
            elif indexOfN > -1 and indexOfX == -1:
                NMap[indexOfN] +=1
            elif indexOfN == -1 and indexOfX > -1:
                NMap[indexOfX] +=1
                
            indexOfN = line.rfind("N")
            indexOfX = line.rfind("X")
            
            if not indexOfN in RNMap:
                RNMap[indexOfN] = 0
                
            if not indexOfX in RNMap:
                RNMap[indexOfX] = 0
                
            if indexOfN == -1 and indexOfX == -1:
                RNMap[-1] += 1
            elif indexOfN > -1 and indexOfX == -1:
                RNMap[indexOfN] +=1
            elif indexOfN == -1 and indexOfX > -1:
                RNMap[indexOfX] +=1
                
        else:
            if trace:
                print(line)
    mf.close()
    return (NMap, RNMap)

def mapToCsv(m, ouputFile):
    of = open(ouputFile, "w")
    keys = sorted(m.keys())
    for k in keys:
        of.write(str(k)+";"+str(m[k]) + "\n")
    of.close()
