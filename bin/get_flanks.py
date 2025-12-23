import sys

###################################################################
## PROGRAM: get_flanks.py                                        ##
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
##    get_flanks.py <fastafile>                                                 ##
##                                                                              ##
## This program reads in a fasta file that contains flanking sequences          ##
## and extracts the flanking sequences with just a single N in between          ##
##################################################################################

f = open(sys.argv[1], 'r')

for line in f:
    l = line.strip()
    if l[0]==">":
        print(l)
    else:
        lsplit = l.split('N')
        for i in range(len(lsplit)):
            if i!=0 and i!=(len(lsplit)-1):
                if lsplit[i] != '':
                    lsplit[i] = 'N'*len(lsplit[i])
        print("N".join(lsplit))
f.close()