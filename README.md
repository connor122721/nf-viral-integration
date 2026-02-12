# Detecting viral integration with: NF-Viral-Integration
- This is a nextflow pipeline for detecting HIV (and or general viral) integration sites using PacBio HiFi sequencing data across a host. Implements the SMRTCap methodology with iterative mapping and multi-reference viral genome support.

## Requirements
- **Nextflow** (≥21.04.0)
- **Container runtime**: Singularity, Apptainer, or Docker

## HPC Usage
For SLURM-based HPC systems:

```bash
# Clone the repository
git clone https://github.com/connor122721/nf-viral-integration.git
cd nf-viral-integration

# Load modules
module load singularity 
module load nextflow

# Run pipeline
nextflow run main.nf \
  --patient_dir 'data/*.fastq.gz' \
  --host_genome /path/to/host.fa \
  --viral_genomes /path/to/hiv_panel.fa \
  --outdir output \
  -profile singularity \
  -resume \
  -bg 
```

Or, to run the pipeline off of GitHub directly use:

```bash
# Load modules
module load singularity 
module load nextflow

# Run pipeline
nextflow run connor122721/nf-viral-integration -latest \ 
  --patient_dir 'data/*.fastq.gz' \
  --host_genome /path/to/host.fa \
  --viral_genomes /path/to/hiv_panel.fa \
  --outdir output \
  -profile singularity \
  -resume \
  -bg
```

## Pipeline Overview

This pipeline detects viral integration sites through:
1. **Iterative mapping** to viral and host genomes
2. **Integration site detection** from chimeric reads
3. **Multi-reference support** for HIV subtypes (A, B, C, D, HIV-2, SIV)

## Input Data

### Required Files
- **PacBio HiFi reads**: FASTQ or BAM format
- **Host genome**: FASTA reference (e.g., hg38)
- **Viral references**: Multi-FASTA with viral genomes

### Example Viral Reference Panel

The pipeline supports multiple HIV references for optimal mapping:

```
HIV-1_subtype_A.fa
HIV-1_subtype_B.fa
HIV-1_subtype_C.fa
HIV-1_subtype_D.fa
HIV-2.fa
SIV.fa
```

## Output Structure
```
output/
├── 01_reference_selection/  # Initial viral alignments
├── 02_iterative_masking/    # Exhaustive viral alignments
├── 03_flank_host_mapping/   # Host genome alignments
└── 04_final_results/        # Detected integrations summary statistics
```

## Key Parameters
```bash
--reads              # Input HiFi reads (FASTQ/BAM)
--host_genome        # Host reference genome (FASTA)
--viral_refs         # Viral reference genome(s) (FASTA)
--outdir             # Output directory [default: ./results]
--min_mapq           # Minimum mapping quality [default: 20]
--min_read_length    # Minimum read length [default: 1000]
```

## Container Profiles
```bash
# Singularity
-profile singularity

# Apptainer
-profile apptainer

# Docker
-profile docker
```

## Pipeline Methodology
Based on SMRTCap protocol for viral integration detection:

1. Filter HiFi reads by quality and length
2. Map to viral reference panel (selects best-matching subtype)
3. Extract unmapped and partially mapped reads
4. Iteratively map to host genome
5. Report integration breakpoints

## Citation
Still in development! If you use this pipeline, please cite:

```bibtex
@software{murray2024nfviralintegration,
  author = {Murray, Connor},
  title = {nf-viral-integration: Nextflow pipeline for viral integration detection},
  year = {2026},
  url = {https://github.com/connor122721/nf-viral-integration}
}
```

## Support
For questions or issues:
- Open an issue on [GitHub](https://github.com/connor122721/nf-viral-integration/issues)
- Email: connor.murray.2@louisville.edu

---

**Develped by Connor Murray, Ph.D.**  
University of Louisville School of Medicine | Dept. of Biochemistry & Molecular Genetics