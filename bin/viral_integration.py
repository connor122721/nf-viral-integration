#!/usr/bin/env python3
"""
Viral Genome Integration Simulator - Lightweight Version
Simulates integration of viral genomes into host genomes
Generates complete modified genome with integrations
"""

import argparse
import random
import sys
from pathlib import Path
from typing import List, Tuple, Optional
import gzip

def load_fasta(fasta_path: Path) -> dict:
    """Load FASTA file into dictionary"""
    sequences = {}
    current_seq = []
    current_name = None
    
    open_func = gzip.open if str(fasta_path).endswith('.gz') else open
    
    with open_func(fasta_path, 'rt') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_name:
                    sequences[current_name] = ''.join(current_seq)
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.upper())
        
        if current_name:
            sequences[current_name] = ''.join(current_seq)
    
    return sequences

def reverse_complement(seq: str) -> str:
    """Return reverse complement of DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))

def generate_integrations(host_genome: dict, 
                         viral_genome: dict,
                         n_integrations: int,
                         target_chroms: Optional[List[str]] = None,
                         seed: Optional[int] = None) -> List[dict]:
    """Generate random viral integration sites"""
    
    if seed:
        random.seed(seed)
    
    if target_chroms is None:
        target_chroms = list(host_genome.keys())
    else:
        target_chroms = [c for c in target_chroms if c in host_genome]
    
    if not target_chroms:
        raise ValueError("No valid target chromosomes found")
    
    viral_name = list(viral_genome.keys())[0]
    viral_len = len(viral_genome[viral_name])
    
    sites = []
    for i in range(n_integrations):
        chrom = random.choice(target_chroms)
        chrom_len = len(host_genome[chrom])
        pos = random.randint(1000, chrom_len - 1000)
        orientation = random.choice(['+', '-'])
        
        sites.append({
            'id': i + 1,
            'chrom': chrom,
            'pos': pos,
            'orientation': orientation,
            'viral_start': 0,
            'viral_end': viral_len
        })
    
    return sites

def integrate_viral_sequence(chrom_seq: str, 
                            viral_seq: str,
                            pos: int,
                            orientation: str,
                            tsd_size: int = 5) -> Tuple[str, int]:
    """
    Integrate viral sequence into chromosome at specified position
    Returns: (modified_sequence, viral_length_inserted)
    """
    
    if orientation == '-':
        viral_seq = reverse_complement(viral_seq)
    
    # Target site duplication
    tsd = chrom_seq[pos:pos + tsd_size]
    
    # Insert: left_part + TSD + viral + TSD + right_part
    modified_seq = chrom_seq[:pos] + tsd + viral_seq + tsd + chrom_seq[pos + tsd_size:]
    
    viral_insert_len = len(tsd) + len(viral_seq) + len(tsd)
    
    return modified_seq, viral_insert_len

def create_integrated_genome(host_genome: dict,
                            viral_genome: dict,
                            sites: List[dict],
                            tsd_size: int = 5) -> Tuple[dict, List[dict]]:
    """
    Create complete genome with all viral integrations
    Returns: (integrated_genome, updated_sites_with_adjusted_positions)
    """
    
    viral_name = list(viral_genome.keys())[0]
    viral_seq = viral_genome[viral_name]
    
    # Group sites by chromosome
    sites_by_chrom = {}
    for site in sites:
        chrom = site['chrom']
        if chrom not in sites_by_chrom:
            sites_by_chrom[chrom] = []
        sites_by_chrom[chrom].append(site)
    
    # Sort sites by position (descending) to integrate from end to start
    # This prevents position shifts from affecting later integrations
    for chrom in sites_by_chrom:
        sites_by_chrom[chrom].sort(key=lambda x: x['pos'], reverse=True)
    
    integrated_genome = {}
    updated_sites = []
    
    # Process each chromosome
    for chrom_name, chrom_seq in host_genome.items():
        if chrom_name in sites_by_chrom:
            # This chromosome has integrations
            modified_seq = chrom_seq
            cumulative_offset = 0
            
            # Process integrations from end to start
            for site in sites_by_chrom[chrom_name]:
                modified_seq, insert_len = integrate_viral_sequence(
                    modified_seq,
                    viral_seq,
                    site['pos'],
                    site['orientation'],
                    tsd_size)
                
                # Track the site with its final position in the modified genome
                updated_site = site.copy()
                updated_site['inserted_length'] = insert_len
                updated_sites.append(updated_site)
            
            integrated_genome[chrom_name] = modified_seq
            print(f"  {chrom_name}: {len(sites_by_chrom[chrom_name])} integration(s), "
                  f"{len(chrom_seq):,} bp -> {len(modified_seq):,} bp")
        else:
            # No integrations, keep original sequence
            integrated_genome[chrom_name] = chrom_seq
    
    # Sort updated sites by original order
    updated_sites.sort(key=lambda x: x['id'])
    
    return integrated_genome, updated_sites

def write_fasta(sequences: dict, output_path: Path):
    """Write sequences to FASTA file"""
    with open(output_path, 'w') as out:
        for name, seq in sequences.items():
            out.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                out.write(f"{seq[i:i+80]}\n")

def main():
    parser = argparse.ArgumentParser(description='Simulate viral genome integration - Generates complete modified genome')
    
    parser.add_argument('--host-genome', type=Path, required=True,
                       help='Host genome FASTA')
    parser.add_argument('--viral-genome', type=Path, required=True,
                       help='Viral genome FASTA')
    parser.add_argument('--n-integrations', type=int, default=10,
                       help='Number of integration sites [default: 10]')
    parser.add_argument('--target-chroms', nargs='+',
                       help='Target chromosomes (default: all)')
    parser.add_argument('--output-fasta', type=Path, required=True,
                       help='Output genome FASTA with integrations')
    parser.add_argument('--output-bed', type=Path, required=True,
                       help='Output integration sites BED')
    parser.add_argument('--seed', type=int,
                       help='Random seed for reproducibility')
    parser.add_argument('--tsd-size', type=int, default=5,
                       help='Target site duplication size [default: 5]')
    
    args = parser.parse_args()
    
    print(f"Loading host genome from {args.host_genome}...")
    host_genome = load_fasta(args.host_genome)
    total_host_size = sum(len(seq) for seq in host_genome.values())
    print(f"  Loaded {len(host_genome)} sequences ({total_host_size:,} bp total)")
    
    print(f"Loading viral genome from {args.viral_genome}...")
    viral_genome = load_fasta(args.viral_genome)
    viral_name = list(viral_genome.keys())[0]
    viral_size = len(viral_genome[viral_name])
    print(f"  Loaded {viral_name} ({viral_size:,} bp)")
    
    print(f"\nGenerating {args.n_integrations} integration sites...")
    sites = generate_integrations(
        host_genome, 
        viral_genome,
        args.n_integrations,
        args.target_chroms,
        args.seed)
    
    print(f"\nIntegrating viral sequences into host genome...")
    integrated_genome, updated_sites = create_integrated_genome(
        host_genome,
        viral_genome,
        sites,
        args.tsd_size)
    
    total_integrated_size = sum(len(seq) for seq in integrated_genome.values())
    size_increase = total_integrated_size - total_host_size
    print(f"\nGenome size: {total_host_size:,} bp -> {total_integrated_size:,} bp "
          f"(+{size_increase:,} bp)")
    
    print(f"\nIntegration sites:")
    for site in updated_sites:
        print(f"  Site {site['id']}: {site['chrom']}:{site['pos']} ({site['orientation']}) "
              f"[+{site['inserted_length']} bp]")
    
    print(f"\nWriting integrated genome to {args.output_fasta}...")
    write_fasta(integrated_genome, args.output_fasta)
    
    print(f"Writing integration sites to {args.output_bed}...")
    with open(args.output_bed, 'w') as bed:
        bed.write("# Integration sites in original genome coordinates\n")
        bed.write("# chrom\tstart\tend\tname\tscore\tstrand\tinserted_length\n")
        for site in updated_sites:
            bed.write(f"{site['chrom']}\t{site['pos']}\t{site['pos']+1}\t"
                     f"integration_{site['id']}\t.\t{site['orientation']}\t"
                     f"{site['inserted_length']}\n")
    
    print("\nDone!")
    print(f"Complete genome with {args.n_integrations} integrations written to {args.output_fasta}")

if __name__ == "__main__":
    main()