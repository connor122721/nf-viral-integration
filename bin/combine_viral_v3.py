#!/usr/bin/env python3
import pysam
import sys
import os
from pathlib import Path
import re
import xlsxwriter
import warnings 
from RTFclass import RTFoutput
from IntegrationClass import IntegrationData
from OUTPUTclass import OUTPUTfile

###################################################################
## PROGRAM: combine_viral_v2.py (Modified for Nextflow)        ##
## AUTHOR:  David Sachs,                                        ##
##          Modified by Eric Rouchka, University of Louisville  ##
##          Modified by Connor S. Murray, University of Louisville ##
## DATE:    01/24/2026                                          ##
##                                                              ##
## USAGE:                                                        ##
##    combine_viral_v2.py <SAMPLE_ID>                           ##
##                                                               ##
## For Nextflow pipeline, expects files in current directory:   ##
##   - sample_id.1.sam, sample_id.2.sam, ... (filtered iteration SAMs) ##
##   - sample_id.flanks.sam                                     ##
##   - sample_id.human.filtered.sam                             ##
##                                                               ##
## The .N.sam files are created by iterative viral mapping:     ##
##   - .1.sam = first iteration (filtered from .initial.1.sam)  ##
##   - .2.sam = second iteration (masked and re-mapped)         ##
##   - etc. (continues until no more viral sequences found)     ##
###################################################################

# Check command line arguments
if len(sys.argv) < 2:
    print("USAGE: combine_viral_v2.py <SAMPLE_ID>")
    print("  The script expects files in current directory:")
    print("    - <SAMPLE_ID>.1.sam (required - first iteration)")
    print("    - <SAMPLE_ID>.2.sam, .3.sam, etc. (optional - additional iterations)")
    print("    - <SAMPLE_ID>.flanks.sam (required)")
    print("    - <SAMPLE_ID>.human.filtered.sam (required)")

## Parse arguments for Nextflow compatibility
sample_id = sys.argv[1]

# For Nextflow, the files are staged in the current working directory
directory = os.getcwd()
prefix = sample_id

print(f"Processing sample: {sample_id}")
print(f"Working directory: {directory}")

# Verify required files exist
required_files = [
    f"{sample_id}.flanks.sam",
    f"{sample_id}.human.filtered.sam"
]

# Check for iteration SAM files (.1.sam, .2.sam, etc.)
# These are the filtered outputs from iterative viral mapping
iteration_sams = []
for i in range(1, 6):  # Check for up to 5 iterations
    iteration_file = f"{sample_id}.{i}.sam"
    if os.path.isfile(iteration_file):
        iteration_sams.append(iteration_file)
        print(f"Found iteration SAM: {iteration_file}")

if not iteration_sams:
    print(f"ERROR: No iteration SAM files found for {sample_id}")
    print(f"Expected at least: {sample_id}.1.sam")

print(f"Found {len(iteration_sams)} iteration SAM file(s)")

# Check other required files
missing_files = []
for req_file in required_files:
    if not os.path.isfile(req_file):
        missing_files.append(req_file)

print("All required files found")

# Initialize the integration data and output

# Initialize IntegrationData with the directory and prefix
# The class should handle finding and opening SAM files
integrations = IntegrationData(directory, prefix)

# Initialize RTF output
RTF = RTFoutput(prefix)
RTF.openrtf(prefix, 0)

# Initialize tab and CSV output files
txtFile = OUTPUTfile(prefix, ".tab", "\t", RTF.columns)
csvFile = OUTPUTfile(prefix, ".csv", ",", RTF.columns)

## Main processing loop
print("Starting main processing loop")
processed_count = 0

try:
    while (not integrations.done):
        integrations.row_data = [""]*(len(RTF.columns))
        integrations.count += 1
            
        # Process each integration site
        integrations.setHIVReads()
        integrations.setFlankNames()
        integrations.setViralInformation(RTF)
        integrations.setHostInformation(RTF)
        integrations.setFlanksInformation(RTF)
        integrations.writeRTF(RTF)  
        integrations.checkHostMappings(RTF)
        integrations.writeFiles(RTF, txtFile, csvFile)
        processed_count += 1
            
        # Progress reporting
        if processed_count % 100 == 0:
            print(f"Processed {processed_count} integration sites...")
except StopIteration:
    pass
finally:
    RTF.closertf()
    txtFile.close()
    csvFile.close()

## Print summary
print(f"\n=== Processing Summary ===")
print(f"Sample: {sample_id}")
print(f"Integration sites processed: {processed_count}")

# List output files
output_files = []
for ext in ['.xlsx', '.tab', '.csv', '.rtf']:
    matches = list(Path(directory).glob(f"{prefix}*{ext}"))
    if matches:
        output_files.extend([f.name for f in matches])

if output_files:
    print(f"\nOutput files created:")
    for f in output_files:
        print(f"  - {f}")
else:
    print(f"WARNING: No output files found")

print(f"\nProcessing complete!")