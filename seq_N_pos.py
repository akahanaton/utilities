#!/usr/bin/python3
from Bio import SeqIO
import os


if __name__ == '__main__':
    import sys
    def usage():
        print('usage:  python3',sys.argv[0],'seq.file', 'nucleotide')
        sys.exit(1)

    #--------------------------------------------------
    # sys.argv = ["../tree.all", "all.longest.fna.aln.fna", "input.tree", "input.seq2", "eriEur"]
    #-------------------------------------------------- 
    argc = len(sys.argv)
    if argc < 3:
        usage()

    seqFile = sys.argv[1]
    nucl = sys.argv[2]

    seqHandle = open(seqFile, "rU")

    lastPos = 0
    for record in SeqIO.parse(seqHandle, "fasta"):
        curSeq=str(record.seq).upper()
        try:
            for i in range(0, len(curSeq)):
                if(curSeq[i] == nucl):
                    if(i == (lastPos +1)):
                        lastPos = i
                        continue
                    else:
                        lastPos = i
                        print(i+1)
        except ValueError:
            print("unkown seq\t", + record.id + "\t" +  curSeq[i])
    seqHandle.close()
