#!/usr/bin/env python3
###################################################################
## PROGRAM: combine_hiv_V2.py                                    ##
## AUTHOR:  David Sachs, Ichan School of Medicine at Mount Sinai ##
##          Modified by Eric Rouchka, University of Louisville   ##
##          Adapted for Nextflow by Connor S. Murray             ##
## DATE:    12/22/2025                                           ##
###################################################################

import sys
import pysam
from collections import defaultdict

def combine_hiv_results(prefix):
    """
    Combine HIV integration results from flank and host alignments.
    
    Args:
        prefix: Sample prefix for finding input files
    """
    # Files to process
    flank_sam = f"{prefix}.flanks.sam"
    host_sam = f"{prefix}.human.filtered.sam"
    output_file = f"{prefix}.integration_sites.txt"
    
    # Data structures
    integration_sites = defaultdict(lambda: {
        'chromosome': None,
        'position': None,
        'strand': None,
        'read_name': None,
        'mapq': 0,
        'confirmed': False
    })
    
    # Parse flank alignments
    print(f"Processing flank alignments from {flank_sam}")
    try:
        with pysam.AlignmentFile(flank_sam, "r") as flank_file:
            for read in flank_file:
                if read.is_unmapped:
                    continue
                
                # Extract integration site information
                read_base = read.query_name.rsplit('_flank', 1)[0] if '_flank' in read.query_name else read.query_name
                
                integration_sites[read_base]['chromosome'] = flank_file.get_reference_name(read.reference_id)
                integration_sites[read_base]['position'] = read.reference_start
                integration_sites[read_base]['strand'] = '-' if read.is_reverse else '+'
                integration_sites[read_base]['read_name'] = read.query_name
                integration_sites[read_base]['mapq'] = read.mapping_quality
    except Exception as e:
        print(f"Warning: Could not process flank alignments: {e}")
    
    # Parse host alignments for confirmation
    print(f"Processing host alignments from {host_sam}")
    try:
        with pysam.AlignmentFile(host_sam, "r") as host_file:
            for read in host_file:
                if read.is_unmapped:
                    continue
                
                read_name = read.query_name
                if read_name in integration_sites:
                    integration_sites[read_name]['confirmed'] = True
    except Exception as e:
        print(f"Warning: Could not process host alignments: {e}")
    
    # Write results
    print(f"\nWriting integration sites to {output_file}")
    with open(output_file, 'w') as out:
        out.write("Read_Name\tChromosome\tPosition\tStrand\tMAPQ\tConfirmed\n")
        
        for read_name, site_info in sorted(integration_sites.items()):
            if site_info['chromosome']:
                out.write(f"{site_info['read_name']}\t"
                         f"{site_info['chromosome']}\t"
                         f"{site_info['position']}\t"
                         f"{site_info['strand']}\t"
                         f"{site_info['mapq']}\t"
                         f"{'Yes' if site_info['confirmed'] else 'No'}\n")
    
    # Summary statistics
    total_sites = len(integration_sites)
    confirmed_sites = sum(1 for site in integration_sites.values() if site['confirmed'])
    
    print(f"\n=== Integration Site Summary ===")
    print(f"Total integration sites detected: {total_sites}")
    print(f"Confirmed by host mapping: {confirmed_sites}")
    print(f"Confirmation rate: {100*confirmed_sites/total_sites:.1f}%" if total_sites > 0 else "N/A")
    
    # Chromosome distribution
    chrom_counts = defaultdict(int)
    for site in integration_sites.values():
        if site['chromosome']:
            chrom_counts[site['chromosome']] += 1
    
    print(f"\n=== Chromosome Distribution ===")
    for chrom, count in sorted(chrom_counts.items(), key=lambda x: -x[1]):
        print(f"{chrom}: {count}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: combine_hiv_V2.py <prefix>")
        sys.exit(1)
    
    combine_hiv_results(sys.argv[1])
