import pysam, sys, re

###################################################################
## PROGRAM: mask.py                                              ##
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
##    mask.py <bamfile>                                                         ##
##                                                                              ##
##################################################################################

BAMFILE = sys.argv[1]
sam = pysam.AlignmentFile(BAMFILE,"r")
count = 0
last_read =  pysam.AlignedSegment()
read_num = 0
coords = [-1, -1]

#-----------------------------------------------------------------------------------------------------------------------------
def process_read(read, last):
    global count, last_read, read_num, coords

    ###############################################################################
    ## Subroutine process_read                                                   ##
    ## IN: read -- current line of sam/bam file                                  ##
    ##---------------------------------------------------------------------------##
    ## This subroutine takes in the current alignment line of the sam/bam file   ##
    ## that has been mapped to the viral (typically HIV/SIV) genome and replaces ##
    ## the mapped portion with N's so it can be iteratively mapped again to the  ##
    ## viral reference until no more hits are found.  At that time, it can be    ##
    ## used to map the flanking regions to the reference host genome             ##
    ###############################################################################

    coords = [read.qstart, read.qend]
    if read.is_reverse:
        coords = [(len(read.seq)-read.qend), len(read.seq)-read.qstart]
    masked = read.get_forward_sequence()[:coords[0]] + 'N'*(coords[1]-coords[0]) + read.get_forward_sequence()[coords[1]:] 
    if read.is_unmapped:
       masked = read.get_forward_sequence()

    #########################################################################
    ## add /0 to forward mapped read; /1 to reverse mapped read descriptor ##
    #########################################################################
    reversed = "0"
    if read.is_reverse:
        reversed = "1"
    print(">" + read.qname + "/" + reversed)
    print(masked)
#-----------------------------------------------------------------------------------------------------------------------------

for read in sam:
    process_read(read, last=False)