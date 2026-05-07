import pysam, sys, re

###################################################################
## PROGRAM: pick_reads.py                                        ##
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
##    pick_reads.py <INPUT SAM FILE 1> <INPUT SAM FILE 2> <OUTPUT SAM FILE>     ##
##                                                                              ##
## Picks reads that had viral integration alignments and need to be mapped      ##
## again to confirm no further alignments are present                           ##
##################################################################################

SAMFILE1 = sys.argv[1]
SAMFILE2 = sys.argv[2]
SAMFILE2_OUT = sys.argv[3]
sam1 = pysam.AlignmentFile(SAMFILE1,"r")
sam2 = pysam.AlignmentFile(SAMFILE2,"r")
sam2_out = pysam.AlignmentFile(SAMFILE2_OUT, "w", template=sam1)
read1num = 0
read2num = 0

read2_tmp = next(sam2)
read2num += 1
qname2 = read2_tmp.qname
ref2 = read2_tmp.reference_name
done = False
while (not done):
    read1 = next(sam1)
    read1num += 1
    qname1 = read1.qname
    ref1 = read1.reference_name
    got_match = False
    while ((qname1.split('/')[:2]==qname2.split('/')[:2]) and (not done)):
        if ((ref2==ref1) or (read2_tmp.is_unmapped)) and (not got_match):
            got_match = True
            read2 = read2_tmp
        try:
            read2_prev = read2_tmp
            read2_tmp = next(sam2)
        except:
            done = True
        read2num += 1
        qname2 = read2_tmp.qname
        ref2 = read2_tmp.reference_name
    if not got_match:
        read2 = read2_prev
        read2.is_unmapped = True
    sam2_out.write(read2)
sam2_out.close()