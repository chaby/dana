# coding: utf8

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
