#!/usr/bin/python
from __future__ import division
from Bio import AlignIO
#--------------------------------------------------
# from Bio.AlignIO import MafIO
#--------------------------------------------------
import sys

def usage():
    print('maf coverage perl alignment')
    print('usage:  python',sys.argv[0],'input.maf species.list')
    sys.exit(1)

argc = len(sys.argv)
if argc < 2:
    usage()
inputMaf = sys.argv[1]
speList = sys.argv[2]

speListH = open(speList,'r')
alllines= speListH.readlines()
speListH.close()

seqList = []
for eachline in alllines:
    spe = eachline.rstrip()
    seqList.append(spe)

for multiple_alignment in AlignIO.parse(inputMaf, "maf"):
    seqCov = dict.fromkeys(seqList,0)
    seqLen = len(multiple_alignment[0].seq)
    for seqrec in multiple_alignment:
        seqCov[seqrec.id.split(".")[0]] = 100 - seqrec.seq.count("-") / seqLen * 100
    print( ("%s\t%s\t%s\t%s" % (multiple_alignment[0].id, multiple_alignment[0].annotations["start"], multiple_alignment[0].annotations["start"] + multiple_alignment[0].annotations["size"] - 1, " ".join(["%s: %0.5s" % (key, ('%('+key+')s') % seqCov) for key in seqList]))))
