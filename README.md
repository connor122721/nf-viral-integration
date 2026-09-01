# nf-viral-integration (SMRTcap: 2026)

- This is a Nextflow pipeline for detecting HIV viral integration sites in PacBio HiFi data
generated on Revio sequencers. Designed for longitudinal patient studies where
the same insertion can be tracked across multiple time-course samples, or standalone identification of inserts.

- This pipeline takes per-sample, demultiplexed Revio output, aligns reads to a
hybrid host + HIV reference, identifies chimeric reads spanning the
host/virus junction, calls integration breakpoints, annotates them against
gene models and RepeatMasker tracks for T2T-CHM13v2.0.

## Quick start

```bash
# Load in singularity (or apptainer)
module load singularity 

# Run pipeline
nextflow run main.nf \
    --samplesheet  samples.csv \
    --host_genome  T2T_refgenome.fasta.gz \
    --annotation   T2T_refgenome.gtf.gz \
    --outdir       results/ \
    -profile       singularity,slurm
```

## Required Input
```
--samples               Metadata sheet csv that organizes input (see below).
   Sample input alternatives:
      --patient_dir           Directory containing patient BAM/FASTQ files.
      --patient_bam           Single patient BAM file.
      --patient_fastq         Single patient FASTQ file (already demultiplexed).
      --multiplexed_fastq     Multiplexed FASTQ(s) to demultiplex with lima.

--host_genome           Path to host genome FASTA.
--viral_genomes         Glob pattern for viral reference panel (e.g., "hiv_refs/*.fa").
```

## Samplesheet (metadata) used for input

`samples.csv` columns:

| column         | required               | notes                                                                                       |
|----------------|------------------------|---------------------------------------------------------------------------------------------|
| sample_id      | yes                    | Unique per row — represents a *final*, post-demultiplex sample                              |
| patient_id     | yes                    | Same value across all samples for one patient                                               |
| timepoint      | yes                    | ISO date or any sortable string (`t0`, `t1`, …)                                             |
| bam            | yes                    | Path to a HiFi BAM. May be multiplexed if `demux=true`.                                     |
| pbi            | optional               | PacBio index; will be created if missing                                                    |
| demux          | optional, default false| `true`/`false`. If `true`, the BAM is run through lima and `barcode` selects this row's output |
| barcode        | required if demux=true | lima output token (e.g. `bc1001--bc1001`). Used to route the per-barcode BAM to this row    |
| barcode_fasta  | required if demux=true | Path to a barcode FASTA understood by lima (or set `params.default_barcode_fasta`)          |

**One row per FINAL sample.** When demultiplexing is needed, write multiple rows
that share the same `bam` path with `demux=true` and unique `barcode` values —
the pipeline groups by parent BAM and runs lima exactly once per multiplexed
input, then dispatches each output to the right downstream sample.

Example mixing both modes:

```csv
sample_id,patient_id,timepoint,bam,demux,barcode,barcode_fasta
PT1_T0,PT1,T0,runs/run1.hifi.bam,true,bc1001--bc1001,refs/Sequel_16_barcodes.fasta
PT1_T1,PT1,T1,runs/run1.hifi.bam,true,bc1002--bc1002,refs/Sequel_16_barcodes.fasta
PT2_T0,PT2,T0,runs/run1.hifi.bam,true,bc1003--bc1003,refs/Sequel_16_barcodes.fasta
PT3_T0,PT3,T0,runs/run2.demuxed.bam,false,,
```

The first three rows trigger one lima invocation on `run1.hifi.bam`; the
fourth row is already demultiplexed and skips lima.

To force the entire pipeline to assume pre-demultiplexed inputs (ignoring
`demux` columns altogether), pass `--skip_demux` on the command line.

## Reference genome download (T2T)

```
# Reference genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz

# Gene annotation file 
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz
```

## Pipeline stages

1. **Optional demultiplexing** (`modules/demux.nf`, `modules/nf-core/lima/`) — runs `lima` once per multiplexed parent BAM
2. **QC** (`bin/qc_mods.nf`) — read length, quality, contamination summaries\
3. **Reference selection** (`bin/select_best_reference.py`) — picks the best
   HIV subtype reference per sample
4. **Alignment & integration calling** (`bin/genomic_processes.nf`,
   `bin/integration_annotation.nf`)
5. **Annotation** (`bin/1.Annotate_Flank_Bam.R`,
   `bin/2.Blast_Viral_Genes.R`, `bin/BMS.insertion.v3.3.ECR_OG.R`)

## Expected Output Layout
This is an example of a simplified output from one sample (`sample1321R`) when mapped to all 9 of the included viral reference genomes (in `${projectDir}/data/hiv_genome_panel`). Outputs were simplified to showcase the general structure. 

- `01_reference_selection`: contains all mapping results to the viral and viral+host genomes across the viral panel.
- `02_phylogenetics`: contains the viral alignments and hypermut3 (hypermutation G --> A statistics) output.
- `final_results`: contains all the annotated viral insertions and BLAST proviral results.
- `converted_fastq` & `lima`: contains all processed fastqs prior to mapping and demultiplexing output.
- `genome_files`: contains the processed GTF used for annotation.
- `pipeline_info`: contains nextflow logs and general metrics of the run so you can track CPU time & memory usage.

After a successful run with `--outdir output/`, the directory tree looks like:

```
output/
├── 01_reference_selection
│   └── sample1321R
│       ├── sample1321R_best_reference.fa
│       ├── sample1321R_best_reference.txt
│       ├── sample1321R_detailed_metrics.txt
│       ├── sample1321R.final.hostflanks.fa
│       ├── sample1321R.final.masked.fa
│       ├── sample1321R.final.unmasked.fa
│       ├── sample1321R_mapping_comparison.txt
│       ├── sample1321R_vs_Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455.allmapped.readnames.txt
│       ├── sample1321R_vs_Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455.dups.readnames.txt
│       ├── sample1321R_vs_Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455.hostflanks.fa
│       ├── sample1321R_vs_Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455.pbmarkdup.log
│       ├── sample1321R_vs_Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455.primary.sam
│       ├── sample1321R_vs_Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455.sorted.bam
│       ├── sample1321R_vs_Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455.sorted_best_reference.sam
│       ├── sample1321R_vs_Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455.sorted.sam
│       ├── sample1321R_vs_Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455.sorted.sam.fa
│       ├── sample1321R_vs_Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455.stats.txt
│       ├── sample1321R_vs_Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455.viralhits.readnames.txt
│       ├── sample1321R_vs_Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455_viralhits.sam
│       ├── sample1321R_vs_Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455.viralreads.fa
│       └── Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455.fasta
├── 02_phylogenetics
│   └── sample1321R
│       ├── sample1321R.final.unmasked.fa.aln
│       ├── sample1321R.final.unmasked.fa.sam.log
│       ├── sample1321R.hypermutargs.csv
│       ├── sample1321R.hypermutpositions.csv
│       ├── sample1321R.hypermutsummary.csv
│       └── viralmsa.log
├── converted_fastq
│   └── sample1321R.fastq.gz
├── failed_best_reference_samples.txt
├── final_results
│   └── sample1321R
│       ├── blast_output
│       │   └── sample1321R.viral.txt_LTR_matches.fa
│       ├── fastas
│       │   ├── sample1321R.final.hostflanks.fa
│       │   ├── sample1321R.final.masked.fa
│       │   └── sample1321R.final.unmasked.fa
│       ├── sample1321R.annotated.csv
│       ├── sample1321R_mapping_comparison.txt
│       └── logs_intermediates
│           ├── CCS_ReadIDs_Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455.txt
│           ├── sample1321R.combined.csv
│           ├── sample1321R.viral.txt
│           ├── sample1321R_vs_Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455.dups.readnames.txt
│           ├── sample1321R_vs_Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455.pbmarkdup.log
│           ├── LTR3_Histogram.png
│           ├── LTR5_Histogram.png
│           ├── LTRBOTH_Histogram.png
│           └── LTRUsage.png
├── genome_files
│   └── host.gtf
├── lima
│   └── sample1321R
│       └── sample1321R.bam
└── pipeline_info
    ├── report.html
    └── timeline.html
```

## Output column dictionaries

### `<sample_id>.annotation.csv` (per-sample integrations/episome identification tables)

| column              | description                                                            |
|---------------------|------------------------------------------------------------------------|
| sample_id           | Sample identifier from the sample sheet                                |
| chrom               | Host chromosome where the integration was called                       |
| position            | 1-based host base coordinate of the breakpoint                         |
| strand              | `+`/`-`/`*`. Strand of the host read at the junction                   |
| viral_orientation   | `+`/`-` orientation of the integrated provirus                         |
| viral_gene          | Annotated viral gene at the integration junction (if BLASTed)          |

## Reproducing repeatmasker step without internet (HPC)

If your nodes can't reach the internet, pre-download the bigBed files once on the
login node and drop them into `params.repeatmasker_cache_dir`. The download
process will see the cached `*.bed.gz` file and skip the process.

```
# Get repeatmasker output for T2T
wget "https://hgdownload.soe.ucsc.edu/gbdb/hs1/t2tRepeatMasker/chm13v2.0_rmsk.align.bb"
genome_label="t2t"

# Convert bigbed to Bedfile!
bigBedToBed chm13v2.0_rmsk.align.bb ${genome_label}.repeatmasker.bed

# Sort + bgzip-friendly gzip 
sort -k1,1 -k2,2n ${genome_label}.repeatmasker.bed | gzip -n > ${genome_label}.repeatmasker.bed.gz
rm -f ${genome_label}.rmsk.bb ${genome_label}.repeatmasker.bed
```

## Demultiplexing details

When `demux=true` rows are present in the samplesheet (and `--skip_demux` is
not set), the pipeline:

1. Parses the samplesheet and tags each row as either pass-through (no demux)
   or demux-required.
2. Groups demux-required rows by their `bam` + `barcode_fasta` pair, so a
   given multiplexed BAM is only fed to lima once even if many final samples
   come from it.
3. Invokes the vendored `nf-core/modules/lima` process (one call per parent
   BAM) with `params.lima_args` (default: `--hifi-preset SYMMETRIC
   --min-score 80`).
4. Routes each per-barcode output BAM back to the correct `sample_id` using a
   barcode→sample_id routing table built from the samplesheet.
5. Concatenates pass-through and demuxed BAMs into a single per-sample
   channel that the rest of the pipeline consumes — so downstream stages do
   not need to know whether a sample was demuxed.

Per-run lima outputs (summary, counts, report) are still published under
`results/<demux_run_id>/lima/` for QC review.

## Citation & Other Resources

### This is the original paper describing the SMRTcap methodology with HIV:
Sadri G, Nadakal ST, Lauer W, Kos J, Singh PK, et al. (2026) Development and validation of HIV SMRTcap for the characterization of HIV-1 reservoirs across tissues and subtypes. PLOS Pathogens 22(1): e1013171. (https://doi.org/10.1371/journal.ppat.1013171)

### This is our pre-print that uses the SMRTcap methodology on lentiviral vectors (LLVs):
Catherine W. Kaiser, Ghazal S. Mehs, Erin M. Elliott, Joanna E. Mroczkowska, Ankita Jain, Michael Ferguson, Alfred L. Garfall, Frederic D. Bushman, Joseph A. Fraietta, Eric C. Rouchka, Melissa L. Smith (2026). LVV SMRTcap reveals extensive proviral variation in lentiviral vector-transduced CAR T cells. bioRxiv. (https://doi.org/10.64898/2026.05.13.724601)

### This is the original bioinformatics pipeline prior to adapation to chimeric read support & nextflow: 
SMRTcap analysis pipeline (2025): (https://github.com/SmithLabLouisville/SMRTCap)