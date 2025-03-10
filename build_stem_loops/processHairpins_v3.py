#!/local/cluster/bin/python
import sys

#################
## SUBROUTINES ##
#################

def readHairpinFile(hairpinFile):
    hairpins = {}
    with open(hairpinFile) as F:
        for line in F:
            rnaID,hpID,hpStart,hpEnd,hpSeq,hpLength,closingBP,isPK = line.strip().split('\t')
            if rnaID not in hairpins:
                hairpins[rnaID] = list()
            hairpins[rnaID].append((hpID,int(hpStart),int(hpEnd),hpSeq))
    return hairpins

def readSegmentFile(segmentFile):
    segments = {}
    with open(segmentFile) as F:
        for line in F:
            rnaID,segmentID,len5p,start5p,end5p,seq5p,len3p,start3p,end3p,seq3p,numBP = line.strip().split('\t')
            if rnaID not in segments:
                segments[rnaID] = list()
            segments[rnaID].append((segmentID,int(start5p),int(end5p),seq5p,int(start3p),int(end3p),seq3p))
    return segments

##########
## MAIN ##
##########

usage = "Usage: " + sys.argv[0] + " <hairpin file> <segment file>"
if len(sys.argv) != 3:
    print(usage)
    sys.exit()

hairpinFile = sys.argv[1]
segmentFile = sys.argv[2]

hairpins = readHairpinFile(hairpinFile)
segments = readSegmentFile(segmentFile)

outputFile = "allHairpins.fasta"
OUTPUT = open(outputFile,'w')

for rnaID in hairpins:
    for hairpin in hairpins[rnaID]:
        hpID,hpStart,hpEnd,hpSeq = hairpin
        for segment in segments[rnaID]:
            segmentID,start5p,end5p,seq5p,start3p,end3p,seq3p = segment
            if end5p + 1 == hpStart:
                # should be a match 
                fullSeq = seq5p + hpSeq + seq3p
                fullID = rnaID + "_" + hpID + "_" + segmentID
                OUTPUT.write(">%s\n%s\n" % (fullID,fullSeq))

OUTPUT.close()
                
