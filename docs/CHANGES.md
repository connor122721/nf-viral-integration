# CHANGES — refactor pass

This document captures every file added, replaced, or marked for retirement
during the optimization pass. Because the bin/ folder could not be uploaded,
some retirement decisions are inferred from file names and your description
and should be confirmed before deletion.

## Added

| path                                           | purpose                                              |
|------------------------------------------------|------------------------------------------------------|
| `bin/assign_clonal_ids.R`                      | Persistent clone IDs across time-course samples      |
| `bin/merge_repeatmasker_annotation.R`          | Joins bedtools output back onto integrations table   |
| `bin/make_circos.R`                            | Per-sample and project circos plots (ggplot2 polar)  |
| `modules/repeatmasker.nf`                      | UCSC bigBed download + bedtools overlap              |
| `modules/post_process.nf`                      | Sub-workflow wiring clonal IDs / circos / HTML       |
| `modules/demux.nf`                             | Optional lima demultiplex sub-workflow               |
| `modules/nf-core/lima/main.nf`                 | Vendored nf-core lima module                         |
| `modules/nf-core/lima/meta.yml`                | Module metadata                                      |
| `examples/samples.example.csv`                 | Example samplesheet showing demux + pass-through     |
| `nextflow.config.additions`                    | New params block (merge into your existing config)   |
| `main_2.nf.snippet`                            | Patch instructions for the main workflow             |
| `README.md`                                    | Full pipeline + output column documentation          |
| `docs/INTEGRATION.md`                          | How to drop these files into your existing repo      |

## Replaced (drop-in)

| path                              | replaces                          | what changed                                                           |
|-----------------------------------|-----------------------------------|------------------------------------------------------------------------|
| `bin/3.Create_Sample_HTML.R`      | original `bin/3.Create_Sample_HTML.R` | External image references, no JS framework, multi-page mode             |
| `bin/4.Create_Project_HTML.R`     | original `bin/4.Create_Project_HTML.R`| Same refactor, plus clonal persistence tables and per-sample links     |

> The new versions keep the same script name. If you want to compare against
> the originals, rename the existing files to `*.legacy.R` before copying the
> new ones into `bin/`.

## Candidates for retirement (please confirm before deleting)

These were inferred from file naming patterns:

1. **`bin/2.findViralGenes.R` vs `bin/2.Blast_Viral_Genes.R`**
   - Two scripts share the `2.` prefix, suggesting one is an older approach.
   - The BLAST-based one (`Blast_Viral_Genes.R`) is almost certainly the
     authoritative version; `findViralGenes.R` is likely an earlier
     k-mer/string-search prototype.
   - **Action:** confirm by diff; if BLAST version is canonical, delete
     `findViralGenes.R` and remove its include from any `.nf` module.

2. **`bin/mask.py` + `bin/unmask.py`**
   - Two scripts that almost certainly share 90% of their code (one applies
     a mask, the other reverses it).
   - **Action:** consolidate into a single `bin/mask.py` with a `--reverse`
     flag. This avoids drift when the mask format is updated in one and not
     the other.

3. **`bin/combine_hiv_V2b.py`**
   - The `V2b` suffix implies earlier `V1`/`V2`/`V2a` versions are in the
     repo's history. Ensure no `.nf` process references an older version
     before deletion. The active version should simply be `combine_hiv.py`
     (rename and update includes).

4. **Legacy HTML embedding utilities**
   - If your old `Create_*_HTML.R` scripts depended on a base64 embedder
     helper script (often named something like `embed_image.R`,
     `make_thumbnails.R`, or similar), it is no longer required after the
     refactor. The new HTMLs reference images by relative path.

5. **Redundant Nextflow include lines in `main_2.nf`**
   - Once the above scripts are removed, prune any `include { ... }` lines
     that now point to deleted files.

## Behavioral changes you should be aware of

- **HTML reports are smaller but split.** Default `params.html_mode = "both"`
  means each run produces *both* a slim single-page report *and* a
  multi-page directory. Set `html_mode = "single"` or `"multi"` to pick one.
- **`repeats/` directory per sample.** New artifact. Existing downstream
  consumers should ignore unfamiliar files unless they enumerate `outdir/`.
- **`integrations.tsv` gains columns.** `clonal_id`, `canonical_pos`,
  `repeat_*_t2t`, `repeat_*_hg38` columns are appended. Anything that
  consumes this TSV by header name continues to work; anything that uses
  positional column indexing must be updated.
- **`projectDir` is referenced** from the Nextflow processes for the
  `Rscript` invocations — this is the standard Nextflow pattern but make
  sure your Singularity containers bind `projectDir` (most do by default).

## Not changed

- Existing alignment / read-picking / integration-calling stages
  (`genomic_processes.nf`, `integration_annotation.nf`, `pick_reads.py`,
  `select_best_reference.py`, `IntegrationClass.py`, `OUTPUTclass.py`,
  `BMS.insertion.v3.3.ECR_OG.R`, `1.Annotate_Flank_Bam.R`,
  `2.Blast_Viral_Genes.R`, `combine_hiv_V2b.py`, `get_flanks.py`,
  `qc_mods.nf`) are untouched. The new sub-workflow consumes their existing
  outputs without altering them.
