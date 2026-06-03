# nf-viral-integration_t2t

A Nextflow pipeline for detecting HIV viral integration sites in PacBio HiFi data
generated on Revio sequencers. Designed for longitudinal patient studies where
the same insertion needs to be tracked across multiple time-course samples.

The pipeline takes per-sample, demultiplexed Revio output, aligns reads to a
hybrid host + HIV reference, identifies chimeric reads spanning the
host/virus junction, calls integration breakpoints, annotates them against
gene models and RepeatMasker tracks for both T2T-CHM13v2.0 and GRCh38, and
produces sample-level and project-level HTML reports including a circos plot
and a clonal-persistence matrix.

## Quick start

```bash
nextflow run main_2.nf \
    --samplesheet  samples.csv \
    --outdir       results/ \
    --report_genome t2t \
    -profile       singularity
```

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

## Pipeline stages

1. **Optional demultiplexing** (`modules/demux.nf`, `modules/nf-core/lima/`)
   ← *new in this refactor* — runs lima once per multiplexed parent BAM
2. **QC** (`bin/qc_mods.nf`) — read length, quality, contamination summaries
3. **Read selection** (`bin/pick_reads.py`) — pulls candidate chimeric reads
4. **Reference selection** (`bin/select_best_reference.py`) — picks the best
   HIV subtype reference per sample
5. **Alignment & integration calling** (`bin/genomic_processes.nf`,
   `bin/integration_annotation.nf`)
6. **Annotation** (`bin/1.Annotate_Flank_Bam.R`,
   `bin/2.Blast_Viral_Genes.R`, `bin/BMS.insertion.v3.3.ECR_OG.R`)
7. **Clonal-ID assignment** (`bin/assign_clonal_ids.R`) ← *new in this refactor*
8. **RepeatMasker overlap** (`modules/repeatmasker.nf`,
   `bin/merge_repeatmasker_annotation.R`) ← *new*
9. **Circos rendering** (`bin/make_circos.R`) ← *new*
10. **HTML reports** (`bin/3.Create_Sample_HTML.R`,
    `bin/4.Create_Project_HTML.R`) — refactored for low file size + multi-page

## Expected outputs

After a successful run with `--outdir results/`, the directory tree looks like:

```
results/
├── <sample_id>/
│   ├── alignment/                      # BAM + index from genomic_processes
│   ├── integrations/
│   │   └── <sample_id>.integrations.tsv      ← raw per-sample calls
│   ├── repeats/
│   │   ├── <sample_id>.t2t.integrations.repeatmasker.tsv
│   │   └── <sample_id>.hg38.integrations.repeatmasker.tsv
│   ├── circos/
│   │   └── <sample_id>.circos.png
│   └── report/
│       ├── <sample_id>.report.html           ← slim single-page
│       └── <sample_id>_report/               ← multi-page version
│           ├── index.html
│           ├── circos.html
│           ├── all.html
│           └── chrom_<chr>.html
├── clonal_tracking/
│   ├── integrations_with_clonal_id.tsv       ← long format, clone-annotated
│   ├── clonal_persistence_wide.tsv           ← clone × timepoint matrix
│   └── clonal_summary.tsv                    ← per-clone summary
└── project/
    ├── circos/
    │   └── project.circos.png
    └── report/
        ├── project.report.html               ← slim single-page
        └── project_report/                   ← multi-page version
            ├── index.html
            ├── circos.html
            ├── samples.html
            ├── persistence.html
            └── summary.html
```

## Output column dictionaries

### `<sample_id>.integrations.tsv` (raw per-sample call table)

| column              | description                                                            |
|---------------------|------------------------------------------------------------------------|
| sample_id           | Sample identifier from the sample sheet                                |
| chrom               | Host chromosome where the integration was called                       |
| position            | 1-based host base coordinate of the breakpoint                         |
| strand              | `+`/`-`/`*`. Strand of the host read at the junction                   |
| viral_orientation   | `+`/`-` orientation of the integrated provirus                         |
| viral_gene          | Annotated viral gene at the integration junction (if BLASTed)          |
| read_support        | Number of HiFi reads supporting the breakpoint                         |
| flank_qual          | Mean Phred Q of the host flank used for the call                       |
| junction_seq        | Soft-clip sequence at the host/virus junction                          |
| integration_class   | Classification (e.g. `gene_body`, `intergenic`, `repeat`)              |

### `integrations_with_clonal_id.tsv` (after `assign_clonal_ids.R`)

Inherits every column above and adds:

| column            | description                                                              |
|-------------------|--------------------------------------------------------------------------|
| clonal_id         | Persistent ID of the form `chrN_<canonical_pos>` (or `..._<strand>`)     |
| canonical_pos     | Median position of all sites in the clone (used for `clonal_id`)         |
| patient_id        | From sample sheet                                                        |
| timepoint         | From sample sheet                                                        |

### `*.integrations.repeatmasker.tsv` (per genome)

Output of `bedtools closest` joined back via `merge_repeatmasker_annotation.R`.
Final columns added to the integrations table:

| column                   | description                                                       |
|--------------------------|-------------------------------------------------------------------|
| repeat_name_<genome>     | RepeatMasker name (e.g. `AluY`, `L1HS`)                           |
| repeat_class_<genome>    | Class parsed from RepeatMasker (`SINE`, `LINE`, `LTR`, ...)       |
| repeat_family_<genome>   | Family parsed from RepeatMasker (`Alu`, `L1`, `ERV1`, ...)        |
| distance_to_repeat_<g>   | Distance in bp to nearest repeat. `0` = sits inside a repeat      |
| in_repeat_<genome>       | Boolean shortcut for `distance_to_repeat_<genome> == 0`           |

`<genome>` is `t2t` or `hg38`. Both pairs of columns coexist on every row so
T2T-only and HG38-only repeats can be compared at a glance.

### `clonal_persistence_wide.tsv`

Wide-format matrix, one row per clone, one column per timepoint. Cell value is
the read support sum across samples sharing that timepoint (or count when
`read_support` is absent).

| column            | description                                            |
|-------------------|--------------------------------------------------------|
| clonal_id         | Persistent clone identifier                            |
| chrom             | Host chromosome                                        |
| canonical_pos     | Median breakpoint position                             |
| `<timepoint_*>`   | Read support at that timepoint, `0` if not detected    |

### `clonal_summary.tsv`

| column            | description                                                              |
|-------------------|--------------------------------------------------------------------------|
| clonal_id         | Persistent clone identifier                                              |
| chrom             | Host chromosome                                                          |
| canonical_pos     | Median breakpoint                                                        |
| n_samples         | Number of distinct samples in which the clone was detected               |
| n_timepoints      | Number of distinct timepoints                                            |
| first_timepoint   | Earliest sortable timepoint at which the clone was seen                  |
| last_timepoint    | Latest                                                                   |
| total_reads       | Sum of `read_support` across all detections                              |
| patient_id        | Comma-separated list of patient IDs (usually one)                        |

Sorted descending by `n_timepoints`, then `n_samples` — i.e. the most
persistent clones first.

## HTML reports — design notes

The reports were rewritten to address the issue raised on GitHub where a
16 MB single-page report crashed browsers on 16 GB-RAM machines:

* **External images.** The circos PNG is referenced via `<img src="…">`
  rather than embedded as base64. PNGs are tens of KB instead of hundreds.
* **No bundled JS frameworks.** Tables use ~40 lines of vanilla JS for
  filtering instead of pulling in DataTables (which inlines ~300 KB).
* **Multi-page mode.** `--html_mode multi` shards the integration table by
  chromosome — each page typically ends up under 200 KB and renders on
  low-RAM laptops.
* **Both modes by default.** `--html_mode both` (the default) generates
  the slim single-page *and* the multi-page version, so users can pick.

Set `params.html_mode = "single"` if you don't want the multi-page tree.

## Parameters added in this refactor

| parameter                  | default                            | meaning                                  |
|----------------------------|------------------------------------|------------------------------------------|
| `clonal_window_bp`         | `10`                               | Window for grouping integrations into clones |
| `clonal_use_strand`        | `false`                            | Require same strand within a clone        |
| `repeatmasker_cache_dir`   | `${projectDir}/refs/repeatmasker`  | Where downloaded BED files live           |
| `repeatmasker_url_t2t`     | UCSC `hs1` rmsk.bb                 | Override if UCSC moves the file           |
| `repeatmasker_url_hg38`    | UCSC `hg38` rmsk.bb                | Override if UCSC moves the file           |
| `skip_repeatmasker`        | `false`                            | Skip download + overlap entirely          |
| `html_mode`                | `both`                             | `single` / `multi` / `both`               |
| `report_genome`            | `hg38`                             | Used by circos: `hg38` or `t2t`           |
| `cytoband_t2t`             | `null`                             | CHM13v2.0 cytoband file for T2T circos    |
| `skip_demux`               | `false`                            | Treat all rows as already demultiplexed   |
| `default_barcode_fasta`    | `null`                             | Fallback barcode FASTA for `demux=true` rows that omit `barcode_fasta` |
| `lima_args`                | `"--hifi-preset SYMMETRIC --min-score 80"` | Pass-through args to lima         |

## Reproducing without internet (HPC)

If your nodes can't reach UCSC, pre-download the bigBed files once on the
login node and drop them into `params.repeatmasker_cache_dir`. The download
process will see the cached `*.bed.gz` file and skip the network call.

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

The lima container is pulled from biocontainers automatically by Singularity;
no SIF rebuild required. If your cluster has no internet egress, mirror
`quay.io/biocontainers/lima:2.9.0--h9ee0642_0` to your local registry and
override the container path in `nextflow.config`.

## Containers / dependencies

The refactor is **dependency-neutral** with respect to your existing SIFs
(`R_Genomics_v1_blast.def`, `Viral_Genomics_v5.def`). Every script uses only
packages already installed in those images. The HTML and circos modules were
specifically rewritten to avoid `htmltools` and `circlize` (neither of which
is in your SIFs):

| script                                | needs                                                 | SIF coverage |
|---------------------------------------|-------------------------------------------------------|--------------|
| `bin/assign_clonal_ids.R`             | `data.table`, `optparse`, `GenomicRanges`             | both         |
| `bin/merge_repeatmasker_annotation.R` | `data.table`, `optparse`                              | both         |
| `bin/make_circos.R`                   | `data.table`, `optparse`, `ggplot2` (via tidyverse)   | both         |
| `bin/3.Create_Sample_HTML.R`          | `data.table`, `optparse` (HTML emitted with base R)   | both         |
| `bin/4.Create_Project_HTML.R`         | `data.table`, `optparse` (HTML emitted with base R)   | both         |
| `modules/repeatmasker.nf`             | `wget`/`curl`, `bedtools`, `bigBedToBed` *or* fallback to `rtracklayer` (already in both SIFs) | both |

* No new R packages, Python packages, or shell tools are introduced.
* The circos plot is rendered with `ggplot2 + coord_polar` (a true circular
  layout) instead of the `circlize` package.
* The HTML reports are emitted with base R `cat()`/`writeLines()` — no
  `htmltools`, no `rmarkdown`, no embedded JS frameworks, just ~40 lines of
  vanilla CSS/JS for filterable tables.
* The `bigBedToBed` binary ships with UCSC kent-tools; if your container
  doesn't bundle it, the process falls back to `rtracklayer::import()` which
  is already in both SIFs.
