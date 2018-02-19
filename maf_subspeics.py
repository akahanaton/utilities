from Bio import AlignIO
#--------------------------------------------------
# from Bio.AlignIO import MafIO
#--------------------------------------------------
import sys

def usage():
    print('get alignment conserved in specific species')
    print('usage:  python',sys.argv[0],'input.maf species.ref species.id1 species.id2')
    sys.exit(1)

argc = len(sys.argv)
if argc < 3:
    usage()
inputMaf = sys.argv[1]
outputMaf = sys.argv[2]
speID = sorted(sys.argv[3:])
outHand = open(outputMaf,"w")
#--------------------------------------------------
# print(speID)
#--------------------------------------------------
for multiple_alignment in AlignIO.parse(inputMaf, "maf"):
    if(len(multiple_alignment) == (len(speID))):
        if(speID == sorted(multiple_alignment.get_all_species_ids()) ):
            #--------------------------------------------------
            # print(sorted(multiple_alignment.get_all_species_ids()))
            #--------------------------------------------------
            AlignIO.write(multiple_alignment, outHand, "maf")

    #--------------------------------------------------
    # for seqrec in multiple_alignment:
    #     print("%20s starts at %s on the %s strand of a sequence %s in length, and runs for %s bp" % \
    #           (seqrec.id,
    #            seqrec.annotations["start"],
    #            seqrec.annotations["strand"],
    #            seqrec.annotations["srcSize"],
    #            seqrec.annotations["size"]))
    #--------------------------------------------------
