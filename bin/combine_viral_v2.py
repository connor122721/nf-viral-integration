#!/usr/bin/env python3

import pysam
import sys
import os
from pathlib import Path
import re
import xlsxwriter
import warnings 

# Import the custom classes
# Note: These need to be available in the same directory or in PYTHONPATH
try:
    from RTFclass import RTFoutput
    from IntegrationClass import IntegrationData
    from OUTPUTclass import OUTPUTfile
except ImportError as e:
    print(f"ERROR: Could not import required classes: {e}", file=sys.stderr)
    print("Make sure RTFclass.py, IntegrationClass.py, and OUTPUTclass.py are in the bin/ directory", file=sys.stderr)
    sys.exit(1)
 
###################################################################
## PROGRAM: combine_hiv_V2b.py (Modified for Nextflow)          ##
## AUTHOR:  David Sachs,                                        ##
##          Modified by Eric Rouchka, University of Louisville  ##
##          Modified by Connor S. Murray, University of Louisville##
## DATE:    01/08/2026                                          ##
###################################################################

###################################################################
## USAGE:                                                        ##
##    combine_hiv_V2b.py <SAMPLE_ID>                            ##
##                                                               ##
## For Nextflow pipeline, expects files in current directory:   ##
##   - sample_id.1.sam, sample_id.2.sam, ... (iteration SAMs)  ##
##   - sample_id.flanks.sam                                     ##
##   - sample_id.human.filtered.sam                             ##
###################################################################

warnings.filterwarnings("error")

# Check command line arguments
if len(sys.argv) < 2:
    print("USAGE: combine_hiv_V2b.py <SAMPLE_ID>", file=sys.stderr)
    print("  The script expects files in current directory:", file=sys.stderr)
    print("    - <SAMPLE_ID>.1.sam (and optionally .2.sam, .3.sam, etc.)", file=sys.stderr)
    print("    - <SAMPLE_ID>.flanks.sam", file=sys.stderr)
    print("    - <SAMPLE_ID>.human.filtered.sam", file=sys.stderr)
    sys.exit(1)

###################################################
## Parse arguments for Nextflow compatibility   ##
###################################################

sample_id = sys.argv[1]

# For Nextflow, the files are staged in the current working directory
# The directory is just the current directory
directory = os.getcwd()
prefix = sample_id

print(f"Processing sample: {sample_id}", file=sys.stderr)
print(f"Working directory: {directory}", file=sys.stderr)

# Verify required files exist
required_files = [
    f"{sample_id}.flanks.sam",
    f"{sample_id}.human.filtered.sam"
]

# Check for at least one iteration SAM file
iteration_sam_found = False
for f in os.listdir(directory):
    if re.match(rf'^{re.escape(sample_id)}\.(\d+)\.sam$', f):
        iteration_sam_found = True
        print(f"Found iteration SAM: {f}", file=sys.stderr)
        break

if not iteration_sam_found:
    print(f"ERROR: No iteration SAM files found for {sample_id}", file=sys.stderr)
    print(f"Expected files like: {sample_id}.1.sam, {sample_id}.2.sam, etc.", file=sys.stderr)
    sys.exit(1)

# Check other required files
missing_files = []
for req_file in required_files:
    if not os.path.isfile(req_file):
        missing_files.append(req_file)
        print(f"ERROR: Required file not found: {req_file}", file=sys.stderr)

if missing_files:
    print(f"ERROR: Missing {len(missing_files)} required file(s)", file=sys.stderr)
    sys.exit(1)

print("All required files found", file=sys.stderr)

###################################################
## Initialize the integration data and output   ##
## files using the custom classes               ##
###################################################

try:
    # Initialize IntegrationData with the directory and prefix
    # The class should handle finding and opening SAM files
    integrations = IntegrationData(directory, prefix)
    print(f"IntegrationData initialized", file=sys.stderr)
except Exception as e:
    print(f"ERROR initializing IntegrationData: {e}", file=sys.stderr)
    import traceback
    traceback.print_exc(file=sys.stderr)
    sys.exit(1)

try:
    # Initialize RTF output
    RTF = RTFoutput(prefix)
    RTF.openrtf(prefix, 0)
    print(f"RTF output initialized", file=sys.stderr)
except Exception as e:
    print(f"ERROR initializing RTF output: {e}", file=sys.stderr)
    import traceback
    traceback.print_exc(file=sys.stderr)
    sys.exit(1)

try:
    # Initialize tab and CSV output files
    txtFile = OUTPUTfile(prefix, ".tab", "\t", RTF.columns)
    csvFile = OUTPUTfile(prefix, ".csv", ",", RTF.columns)
    print(f"Output files initialized", file=sys.stderr)
except Exception as e:
    print(f"ERROR initializing output files: {e}", file=sys.stderr)
    import traceback
    traceback.print_exc(file=sys.stderr)
    sys.exit(1)

###################################################
## Main processing loop                         ##
###################################################

print("Starting main processing loop", file=sys.stderr)
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
            print(f"Processed {processed_count} integration sites...", file=sys.stderr)
            
except StopIteration:
    print(f"Finished processing (StopIteration)", file=sys.stderr)
except Exception as e:
    print(f"ERROR during processing: {e}", file=sys.stderr)
    import traceback
    traceback.print_exc(file=sys.stderr)
    # Continue to close files properly
finally:
    # Close all output files
    try:
        RTF.closertf()
        txtFile.close()
        csvFile.close()
        print(f"All output files closed", file=sys.stderr)
    except Exception as e:
        print(f"ERROR closing files: {e}", file=sys.stderr)

###################################################
## Print summary                                ##
###################################################

print(f"\n=== Processing Summary ===", file=sys.stderr)
print(f"Sample: {sample_id}", file=sys.stderr)
print(f"Integration sites processed: {processed_count}", file=sys.stderr)

# List output files
output_files = []
for ext in ['.xlsx', '.tab', '.csv', '.rtf']:
    matches = list(Path(directory).glob(f"{prefix}*{ext}"))
    if matches:
        output_files.extend([f.name for f in matches])

if output_files:
    print(f"\nOutput files created:", file=sys.stderr)
    for f in output_files:
        print(f"  - {f}", file=sys.stderr)
else:
    print(f"WARNING: No output files found", file=sys.stderr)

print(f"\nProcessing complete!", file=sys.stderr)
sys.exit(0)

# Note: The original script had this at the end, which appears unreachable
# integrations.testZMW()
