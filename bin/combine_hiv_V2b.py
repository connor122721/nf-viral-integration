import pysam, sys
from os import listdir
from os.path import isfile, join
import re
import xlsxwriter
import warnings 
from RTFclass import RTFoutput
from IntegrationClass import IntegrationData
from OUTPUTclass import OUTPUTfile
 
#################################################################
## PROGRAM: combine_hiv_V2b.py                                 ##
## AUTHOR:  David Sachs,                                       ##
##          Modified by Eric Rouchka, University of Louisville ##
##          Copyright, University of Louisville                ##
## DATE:    04/23/2025                                         ##
## LAST MODIFIIED: 05/01/2025                                  ##
## GITHUB:                                                     ##
## LICENSE:                                                    ##
#################################################################

#################################################################
## USAGE:                                                      ##
##    combine_hiv_V2.py <PREFIX>                               ##
##                                                             ##
## NEED TO REFORMAT AND CREATE SUBROUTINES                     ##
#################################################################

#================================================================
warnings.filterwarnings("error")

###################################################
## Get the directory information from the prefix ##
## command line argument                         ##
###################################################

directory_prefix = sys.argv[1]
directory = "/".join(directory_prefix.split('/')[:-1])
prefix = directory_prefix.split('/')[-1]

integrations = IntegrationData(directory, prefix)
RTF = RTFoutput(prefix)
RTF.openrtf(prefix, 0)
txtFile = OUTPUTfile(prefix, ".tab", "\t", RTF.columns)
csvFile = OUTPUTfile(prefix, ".csv", ",", RTF.columns)

while (not integrations.done):
   integrations.row_data = [""]*(len(RTF.columns))
   integrations.count += 1
   integrations.setHIVReads()
   integrations.setFlankNames()
   integrations.setViralInformation(RTF)
   integrations.setHostInformation(RTF)
   integrations.setFlanksInformation(RTF)
   integrations.writeRTF(RTF)  
   integrations.checkHostMappings(RTF)
   integrations.writeFiles(RTF, txtFile, csvFile)

RTF.closertf()
txtFile.close()
csvFile.close()
sys.exit()
#integrations.testZMW()
