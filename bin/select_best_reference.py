#!/usr/bin/env python3
"""
Select best viral reference genome based on detailed alignment metrics
"""

import re
import sys
import argparse
from pathlib import Path
import shutil
from collections import defaultdict

def parse_sam_for_metrics(sam_file):
    """Extract detailed alignment metrics from SAM file"""
    metrics = {
        'mapped_reads': 0,
        'total_edit_distance': 0,
        'total_mismatches': 0,
        'total_indels': 0,
        'total_aligned_length': 0,
        'total_mapq': 0,
        'perfect_matches': 0,
        'high_quality_matches': 0  # MAPQ >= 40
    }
    
    with open(sam_file) as f:
        for line in f:
            if line.startswith('@'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 11:
                continue
            
            flag = int(fields[1])
            # Skip unmapped reads
            if flag & 4:
                continue
            
            mapq = int(fields[4])
            cigar = fields[5]
            
            metrics['mapped_reads'] += 1
            metrics['total_mapq'] += mapq
            
            if mapq >= 40:
                metrics['high_quality_matches'] += 1
            
            # Parse optional fields for NM (edit distance)
            nm = None
            for field in fields[11:]:
                if field.startswith('NM:i:'):
                    nm = int(field.split(':')[2])
                    break
            
            if nm is not None:
                metrics['total_edit_distance'] += nm
                if nm == 0:
                    metrics['perfect_matches'] += 1
            
            # Parse CIGAR for indels
            cigar_ops = re.findall(r'(\d+)([MIDNSHPX=])', cigar)
            aligned_length = 0
            for length, op in cigar_ops:
                length = int(length)
                if op in ['M', '=', 'X']:
                    aligned_length += length
                elif op == 'I' or op == 'D':
                    metrics['total_indels'] += length
            
            metrics['total_aligned_length'] += aligned_length
    
    return metrics

def extract_ref_name(filename):
    """Extract reference name from filename pattern: reads_vs_refname.sam"""
    match = re.search(r'_vs_(.+)\.sam$', filename)
    if match:
        return match.group(1)
    
    # Try alternative pattern
    match = re.search(r'_vs_(.+)\.stats\.txt$', filename)
    if match:
        return match.group(1)
    
    return None

def build_genome_lookup(viral_genomes):
    """Build comprehensive lookup dictionary for viral genome files"""
    genome_paths = {}
    
    for genome in viral_genomes:
        genome_path = Path(genome).resolve()  # Get absolute path
        
        if not genome_path.exists():
            print(f"Warning: Genome file does not exist: {genome_path}", file=sys.stderr)
            continue
        
        # Store with full filename (without directory)
        genome_paths[genome_path.name] = str(genome_path)
        
        # Store with stem only (no extension)
        genome_paths[genome_path.stem] = str(genome_path)
        
        # Handle double extensions like .fa.gz, .fasta.gz
        if genome_path.suffix in ['.gz', '.bz2']:
            base_stem = Path(genome_path.stem).stem
            genome_paths[base_stem] = str(genome_path)
        
        # Also store without common FASTA extensions
        base_name = re.sub(r'\.(fa|fasta|fna|ffn)(\.gz|\.bz2)?$', '', genome_path.name, flags=re.IGNORECASE)
        genome_paths[base_name] = str(genome_path)
    
    return genome_paths

def find_best_match(ref_name, genome_paths):
    """Find the best matching genome file for a reference name"""
    
    # Direct exact match
    if ref_name in genome_paths:
        return genome_paths[ref_name]
    
    # Case-insensitive match
    ref_lower = ref_name.lower()
    for key, path in genome_paths.items():
        if key.lower() == ref_lower:
            return path
    
    # Partial match
    for key, path in genome_paths.items():
        if ref_name in key or key in ref_name:
            return path
    
    return None

def main():
    parser = argparse.ArgumentParser(
        description='Select best viral reference based on alignment metrics'
    )
    parser.add_argument('--sam-files', nargs='+', required=True,
                       help='SAM alignment files')
    parser.add_argument('--stats-files', nargs='+', required=True,
                       help='Stats files from mapping')
    parser.add_argument('--viral-genomes', nargs='+', required=True,
                       help='Viral genome FASTA files')
    parser.add_argument('--output-best', default='best_reference.txt',
                       help='Output file for best reference name')
    parser.add_argument('--output-fa', default='best_reference.fa',
                       help='Output FASTA for best reference')
    parser.add_argument('--output-comparison', default='mapping_comparison.txt',
                       help='Output comparison table')
    parser.add_argument('--output-detailed', default='detailed_metrics.txt',
                       help='Output detailed metrics')
    parser.add_argument('--debug', action='store_true',
                       help='Print debug information')
    
    args = parser.parse_args()
    
    # Build genome path lookup
    genome_paths = build_genome_lookup(args.viral_genomes)
    
    if args.debug:
        print(f"DEBUG: Found {len(genome_paths)} genome entries:", file=sys.stderr)
        for key in sorted(genome_paths.keys())[:10]:  # Show first 10
            print(f"  {key} -> {genome_paths[key]}", file=sys.stderr)
        print("", file=sys.stderr)
    
    results = {}
    
    # Process each SAM file
    for sam_file in args.sam_files:
        sam_path = Path(sam_file)
        ref_name = extract_ref_name(sam_path.name)
        
        if not ref_name:
            print(f"Warning: Could not extract reference name from {sam_path.name}", 
                  file=sys.stderr)
            continue
        
        # Parse detailed metrics from SAM
        sam_metrics = parse_sam_for_metrics(sam_file)
        
        if sam_metrics['mapped_reads'] == 0:
            print(f"Warning: No mapped reads for reference {ref_name}", file=sys.stderr)
            continue
        
        # Calculate averages
        avg_edit_distance = sam_metrics['total_edit_distance'] / sam_metrics['mapped_reads']
        avg_mapq = sam_metrics['total_mapq'] / sam_metrics['mapped_reads']
        perfect_match_rate = sam_metrics['perfect_matches'] / sam_metrics['mapped_reads']
        high_quality_rate = sam_metrics['high_quality_matches'] / sam_metrics['mapped_reads']
        
        # Multi-criteria scoring system
        score = (
            sam_metrics['mapped_reads'] * 1000 +          # Primary: number of mapped reads
            (100 - avg_edit_distance) * 100 +             # Secondary: sequence similarity
            avg_mapq * 10 +                                # Tertiary: mapping quality
            perfect_match_rate * 500 +                     # Bonus: perfect matches
            high_quality_rate * 200                        # Bonus: high quality alignments
        )
        
        results[ref_name] = {
            'mapped_reads': sam_metrics['mapped_reads'],
            'avg_edit_distance': avg_edit_distance,
            'avg_mapq': avg_mapq,
            'perfect_matches': sam_metrics['perfect_matches'],
            'perfect_match_rate': perfect_match_rate * 100,
            'high_quality_matches': sam_metrics['high_quality_matches'],
            'high_quality_rate': high_quality_rate * 100,
            'total_indels': sam_metrics['total_indels'],
            'score': score
        }
    
    if not results:
        print("ERROR: No valid mapping results found", file=sys.stderr)
        with open(args.output_best, 'w') as f:
            f.write("NONE\n")
        with open(args.output_fa, 'w') as f:
            f.write("# No suitable reference found\n")
        sys.exit(1)
    
    # Write comparison table
    with open(args.output_comparison, 'w') as f:
        f.write("Reference\tMapped_Reads\tAvg_Edit_Dist\tAvg_MAPQ\tPerfect_Matches\tScore\n")
        for ref, stats in sorted(results.items(), key=lambda x: x[1]['score'], reverse=True):
            f.write(f"{ref}\t{stats['mapped_reads']}\t"
                   f"{stats['avg_edit_distance']:.2f}\t{stats['avg_mapq']:.2f}\t"
                   f"{stats['perfect_matches']}\t{stats['score']:.2f}\n")
    
    # Write detailed metrics
    with open(args.output_detailed, 'w') as f:
        f.write("=== DETAILED REFERENCE COMPARISON ===\n\n")
        
        for ref, stats in sorted(results.items(), key=lambda x: x[1]['score'], reverse=True):
            f.write(f"Reference: {ref}\n")
            f.write(f"  Mapped reads:           {stats['mapped_reads']}\n")
            f.write(f"  Avg edit distance:      {stats['avg_edit_distance']:.3f}\n")
            f.write(f"  Avg MAPQ:               {stats['avg_mapq']:.2f}\n")
            f.write(f"  Perfect matches:        {stats['perfect_matches']} ({stats['perfect_match_rate']:.1f}%)\n")
            f.write(f"  High quality (Q>=40):   {stats['high_quality_matches']} ({stats['high_quality_rate']:.1f}%)\n")
            f.write(f"  Total indels:           {stats['total_indels']}\n")
            f.write(f"  Composite score:        {stats['score']:.2f}\n")
            f.write(f"\n")
    
    # Select best reference
    sorted_results = sorted(results.items(), key=lambda x: x[1]['score'], reverse=True)
    best = sorted_results[0]
    best_ref_name = best[0]
    best_score = best[1]['score']
    
    # Check for ties (within 1% of best score)
    ties = [r for r in sorted_results if abs(r[1]['score'] - best_score) / best_score < 0.01]
    
    # Write best reference name
    with open(args.output_best, 'w') as f:
        f.write(best_ref_name + '\n')
        if len(ties) > 1:
            f.write(f"# Note: {len(ties)} references had very similar scores:\n")
            for ref, stats in ties:
                f.write(f"#   {ref}: {stats['score']:.2f}\n")
    
    # Find and copy the best reference file
    src_path = find_best_match(best_ref_name, genome_paths)
    
    if src_path and Path(src_path).exists():
        if args.debug:
            print(f"DEBUG: Copying {src_path} to {args.output_fa}", file=sys.stderr)
        shutil.copy(src_path, args.output_fa)
        copied = True
    else:
        # Create placeholder with diagnostic information
        with open(args.output_fa, 'w') as f:
            f.write(f"# Best reference: {best_ref_name}\n")
            f.write(f"# ERROR: Original file not found\n")
            f.write(f"# Available genome keys:\n")
            for key in sorted(genome_paths.keys())[:20]:
                f.write(f"#   {key}\n")
        print(f"ERROR: Could not locate reference file for '{best_ref_name}'", file=sys.stderr)
        print(f"Available references: {', '.join(sorted(list(genome_paths.keys())[:10]))}", file=sys.stderr)
        copied = False

    # Print summary
    print(f"\n{'='*60}")
    print(f"BEST REFERENCE SELECTED: {best_ref_name}")
    print(f"{'='*60}")
    print(f"  Mapped reads:        {best[1]['mapped_reads']}")
    print(f"  Avg edit distance:   {best[1]['avg_edit_distance']:.3f}")
    print(f"  Avg MAPQ:            {best[1]['avg_mapq']:.2f}")
    print(f"  Perfect matches:     {best[1]['perfect_matches']} ({best[1]['perfect_match_rate']:.1f}%)")
    print(f"  Composite score:     {best[1]['score']:.2f}")
    
    if copied:
        print(f"  Reference file:      {src_path}")
    else:
        print(f"  Reference file:      NOT FOUND")
    
    if len(ties) > 1:
        print(f"\n  WARNING: {len(ties)} references had very similar scores (within 1%):")
        for ref, stats in ties[:5]:  # Show top 5
            print(f"    - {ref}: score={stats['score']:.2f}, edit_dist={stats['avg_edit_distance']:.3f}")
    print(f"{'='*60}\n")
    
    if not copied:
        sys.exit(1)


if __name__ == "__main__":
    main()