import pysam, sys, re

###################################################################
## PROGRAM: unmask.py                                            ##
## AUTHOR:  David Sachs, Ichan School of Medicine at Mount Sinai ##
##          Modified by Eric Rouchka, University of Louisville   ##
##          Copyright, University of Louisville                  ##
## DATE:    4/23/2025                                            ##
## LAST MODIFIIED: 4/23/2025                                     ##
## GITHUB:                                                       ##
## LICENSE:                                                      ##
###################################################################

##################################################################################
## USAGE:                                                                       ##
##    unmask.py <bamfile> <fastafile>                                           ##
##                                                                              ##
## This program reads in a sam file line by line and adds back in the masked    ##
## portions of the reads, which could either be the inserted viral sequence or  ##
## the flanking host, which is determined by the presence of an "N" in the      ##
## alignment that is a non-N in the actual fasta sequence                       ##
##################################################################################


BAMFILE = sys.argv[1]
sam = pysam.AlignmentFile(BAMFILE,"r")
count = 0
last_read =  pysam.AlignedSegment()
read_num = 0
coords = [-1, -1]
maskfile = open(sys.argv[2], 'r')

#----------------------------------------------------------------------------------
def process_read(read, last):
    global count, last_read, read_num, coords

    ############################################################### 
    ## Function process_read                                     ##
    ## INPUT: read: current sam alignment read                   ##
    ##        last: flag indicating if this is the last sequence ##
    ##-----------------------------------------------------------##
    ## Processes the SAM file and fasta file to unmask the       ##
    ## sequences n the fasta file to produce a new unmasked fa   ##
    ###############################################################

    mask = maskfile.readline()
    mask = maskfile.readline()
    fullseq = read.get_forward_sequence()
    unmasked = ""
    for i in range(len(mask)):
        if mask[i]=="N":
            unmasked = unmasked + fullseq[i]
        else:
            unmasked = unmasked + "N"
    reversed = "0"
    if read.is_reverse:
        reversed = "1"
    print(">" + read.qname + "/" + reversed)
    print(unmasked)
#----------------------------------------------------------------------------------

for read in sam:
    process_read(read, last=False)
maskfile.close()