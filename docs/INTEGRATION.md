# INTEGRATION.md — how to drop these files into your repo

## File map

Copy from this output bundle into your existing `nf-viral-integration_t2t/`:

```
bin/assign_clonal_ids.R                  -> bin/assign_clonal_ids.R
bin/merge_repeatmasker_annotation.R      -> bin/merge_repeatmasker_annotation.R
bin/make_circos.R                        -> bin/make_circos.R
bin/3.Create_Sample_HTML.R               -> bin/3.Create_Sample_HTML.R   (REPLACE)
bin/4.Create_Project_HTML.R              -> bin/4.Create_Project_HTML.R  (REPLACE)
modules/repeatmasker.nf                  -> modules/repeatmasker.nf      (NEW dir)
modules/post_process.nf                  -> modules/post_process.nf
modules/demux.nf                         -> modules/demux.nf             (NEW)
modules/nf-core/lima/main.nf             -> modules/nf-core/lima/main.nf (NEW; vendored)
modules/nf-core/lima/meta.yml            -> modules/nf-core/lima/meta.yml
examples/samples.example.csv             -> examples/samples.example.csv
README.md                                -> README.md                    (REPLACE)
docs/CHANGES.md                          -> docs/CHANGES.md
```

Then merge `nextflow.config.additions` into your existing `nextflow.config`
(append the new params and the `process { withName: ... }` blocks).

Patch `main_2.nf` per `main_2.nf.snippet`:

1. Add `include { POST_PROCESS } from './modules/post_process'` at the top.
2. After your existing per-sample integrations channel is built, call:
   ```nextflow
   POST_PROCESS(integrations_per_sample_ch,
                Channel.fromPath(params.samplesheet))
   ```

## Smoke test (no real data needed)

Validate the R scripts standalone before re-running the full pipeline:

```bash
# 1. Build a tiny synthetic integrations TSV
cat > /tmp/intg.tsv <<'EOF'
sample_id	chrom	position	strand	read_support	timepoint	patient_id	viral_orientation
s1_t0	chr7	1234567	+	12	t0	p1	+
s1_t1	chr7	1234571	+	8	t1	p1	+
s1_t2	chr7	1234565	+	14	t2	p1	+
s1_t0	chr2	98765432	-	5	t0	p1	-
s2_t0	chr12	5556677	+	3	t0	p2	+
EOF

# 2. Sample sheet
cat > /tmp/ss.tsv <<'EOF'
sample_id	patient_id	timepoint
s1_t0	p1	t0
s1_t1	p1	t1
s1_t2	p1	t2
s2_t0	p2	t0
EOF

# 3. Clonal IDs  (chr7 sites within 10bp should collapse to one clone)
Rscript bin/assign_clonal_ids.R \
    --integrations /tmp/intg.tsv \
    --samplesheet  /tmp/ss.tsv \
    --window 10 \
    --out_long    /tmp/intg.clonal.tsv \
    --out_wide    /tmp/intg.wide.tsv \
    --out_summary /tmp/intg.summary.tsv

cut -f1,2,3,4,5 /tmp/intg.clonal.tsv | head     # expect chr7_1234568 to repeat 3x
cat /tmp/intg.summary.tsv                       # expect that clone with n_timepoints=3

# 4. Circos (sample mode) — sanity check that PNG comes out
Rscript bin/make_circos.R \
    --integrations /tmp/intg.clonal.tsv \
    --mode sample --genome hg38 \
    --label "smoke test" --out /tmp/sample.circos.png
file /tmp/sample.circos.png                     # expect PNG image data

# 5. Sample HTML
mkdir -p /tmp/report
Rscript bin/3.Create_Sample_HTML.R \
    --integrations /tmp/intg.clonal.tsv \
    --circos_png  /tmp/sample.circos.png \
    --sample_id   s1_t0 \
    --html_mode   both \
    --outdir      /tmp/report
ls -lh /tmp/report                              # expect <2 MB single + multi/ dir
```

If those four steps succeed, the Nextflow integration will work because the
sub-workflow simply orchestrates the same Rscript calls with paths from your
channels.

## Things to confirm in your environment

Confirmed against your `R_Genomics_v1_blast.def` and `Viral_Genomics_v5.def`:

- [x] **R packages used:** `data.table`, `optparse`, `ggplot2` (via tidyverse),
      `GenomicRanges`, `rtracklayer` — all present in both SIFs.
- [x] **No `circlize`** required — `make_circos.R` was rewritten with
      `ggplot2 + coord_polar` since `circlize` is not in either SIF.
- [x] **No `htmltools`** required — the HTML scripts emit raw HTML with base R.
- [ ] **`bigBedToBed`** binary (UCSC kent tools) in your bedtools/genomics
      container? Most do; if not, the R `rtracklayer` fallback covers you.
      Check with `apptainer exec your.sif which bigBedToBed`.
- [ ] Cluster nodes can reach `hgdownload.soe.ucsc.edu`, *or* you have
      pre-staged the bigBed/bed.gz files into `params.repeatmasker_cache_dir`.
- [ ] `params.samplesheet` is exposed in your existing config (most pipelines
      already have this).

### Which SIF should run which step?

The post-processing scripts only use `data.table`, `optparse`, `ggplot2`,
`GenomicRanges`, and `rtracklayer`, all of which are in both
`R_Genomics_v1_blast.def` (rocker-based) and `Viral_Genomics_v5.def`
(conda/mamba-based). Pick whichever you already use for your other R-based
processes — both will work. In your `nextflow.config`, set:

```nextflow
process {
    withName: 'ASSIGN_CLONAL_IDS|CIRCOS_SAMPLE|CIRCOS_PROJECT|SAMPLE_HTML|PROJECT_HTML' {
        container = 'path/to/viral_int_R4.sif'   // or r-genomics.sif
    }
    withName: 'REPEATMASKER_DOWNLOAD|REPEATMASKER_OVERLAP' {
        container = 'path/to/your_bedtools.sif'  // whichever has bedtools + bigBedToBed
    }
}
```

## Rollback

The refactor never deletes anything. To revert:

1. Remove the include line in `main_2.nf`.
2. Restore your originals of `3.Create_Sample_HTML.R` and
   `4.Create_Project_HTML.R` from git.
3. Delete `modules/repeatmasker.nf`, `modules/post_process.nf`, the new
   `bin/assign_clonal_ids.R`, `bin/merge_repeatmasker_annotation.R`,
   `bin/make_circos.R`, and the appended `params { ... }` block.

Nothing else is mutated.
