#!/usr/bin/env python3

from Bio.Seq import Seq
import os

path="/SAN/Susanas_den/EimeriaMicrobiome/amplicon/"
path2="/SAN/Susanas_den/EimeriaMicrobiome/primer/"
#os.mkdir(path)
#os.mkdir(path2)

with open("/SAN/Victors_playground/Eimeria_microbiome/Multimarker/primer.file.multi.csv", "r") as primerlist:
    next(primerlist)
    for line in primerlist:
        line = line.rstrip()
        primerF,seqF,primerR,seqR=line.split(",")
        #os.mkdir(path + primerF + primerR)
        seqF=Seq(seqF)
        seqR=Seq(seqR)
        RRC=seqR.reverse_complement()
        FRC=seqF.reverse_complement()
        a=("-a %s...%s" % (seqF,RRC))
        A=("-A %s...%s" % (seqR, FRC))
        f= open(path2+primerF+primerR, "w")
        f.write(a + " " + A)
        f.close
        
#save each line is a text file in a  new folder
