// modules/demux.nf
// ============================================================================
// Optional lima-based demultiplexing of multiplexed PacBio HiFi reads.
//
// Handles BOTH multiplexed BAMs and multiplexed FASTQs. The nf-core lima
// module picks its output extension from the input extension, so:
//
//     run.bam       -> <prefix>.<barcode>.bam
//     run.fastq.gz  -> <prefix>.<barcode>.fastq.gz
//
// lima cannot write BAM from FASTQ input (see the in/out compatibility matrix
// at lima.how/get-started), so FASTQ inputs stay FASTQ all the way through and
// are emitted on a separate channel from the BAMs.
//
// Samplesheet schema (one row per FINAL sample after demultiplexing):
//
//   sample_id, patient_id, timepoint, bam, fastq, demux, barcode, barcode_fasta
//   PT1_T0,    PT1,        T0,        run1.bam, ,      true,  bc1001--bc1001, refs/bc.fa
//   PT1_T1,    PT1,        T1,        run1.bam, ,      true,  bc1002--bc1002, refs/bc.fa
//   PT2_T0,    PT2,        T0,        run2.bam, ,      false, ,
//   PT3_T0,    PT3,        T0,        ,    pool.fastq.gz, true,  bc1003--bc1003, refs/bc.fa
//   PT4_T0,    PT4,        T0,        ,    solo.fastq.gz, false, ,
//
//   - Rows with demux == 'true'  share a reads file and run lima ONCE on it;
//     each row's `barcode` selects the corresponding lima output.
//   - Rows with demux == 'false' (or missing/empty) bypass lima entirely.
//   - A row may set `barcode` to '*' to mean "emit every barcode lima finds",
//     naming samples <parent>.<barcode>. Used by --multiplexed_fastq.
//   - If a row populates both `bam` and `fastq`, the BAM wins.
//
// Public lima module is vendored at modules/nf-core/lima/main.nf.
// ============================================================================

include { LIMA } from './nf-core/lima/main'

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

def dmxHasVal(row, col) {
    row.containsKey(col) && row[col] != null && row[col].toString().trim() != ''
}

// Which file does this row actually read from? BAM wins if both are given.
def readsForRow(row) {
    if (dmxHasVal(row, 'bam') && dmxHasVal(row, 'fastq')) {
        log.warn "[DEMUX] row '${row.sample_id}' populates both 'bam' and 'fastq'; using the BAM."
    }
    if (dmxHasVal(row, 'bam'))   return file(row.bam,   checkIfExists: true)
    if (dmxHasVal(row, 'fastq')) return file(row.fastq, checkIfExists: true)
    throw new IllegalArgumentException(
        "[DEMUX] row '${row.sample_id}' has neither a 'bam' nor a 'fastq' value.")
}

def isFastqFile(f) {
    f.getName() ==~ /(?i).*\.(fastq|fq)(\.gz)?$/
}

// Strip the read-file extension, then take the trailing dot-field, which is
// the barcode token lima's --split-named writes:
//   demux_run1.bc1001--bc1001.fastq.gz -> bc1001--bc1001
def barcodeFromName(name) {
    def stem = name.replaceAll(/(?i)\.(bam|fastq|fq|fasta|fa)(\.gz)?$/, '')
    stem.contains('.') ? stem.substring(stem.lastIndexOf('.') + 1) : stem
}

// Turn one lima output group into per-final-sample tuples. Shared by the BAM
// and FASTQ paths so routing behaves identically for both.
def routeLimaOutputs(meta, outfiles) {
    def flist = (outfiles instanceof List) ? outfiles : [outfiles]

    flist.collect { f ->
        def barcode = barcodeFromName(f.getName())

        // Auto-route mode: no explicit barcode->sample table, so the barcode
        // itself names the sample.
        if (meta.auto_route) {
            def auto_id = "${meta.parent_id}.${barcode}".toString()
            return tuple([
                id:         auto_id,
                sample_id:  auto_id,
                patient:    null,
                timepoint:  null,
                demuxed:    true,
                barcode:    barcode,
                parent:     meta.parent,
                parent_bam: meta.parent
            ], f)
        }

        def sample_id = meta.routing[barcode]
        if (!sample_id) {
            log.warn "[DEMUX] lima produced output for unmapped barcode '${barcode}' " +
                     "from parent ${meta.parent} — skipping. " +
                     "(routing keys: ${meta.routing.keySet()})"
            return null
        }

        def smeta = meta.samples_meta[sample_id]
        tuple([
            id:         sample_id,
            sample_id:  sample_id,
            patient:    smeta?.patient_id,
            timepoint:  smeta?.timepoint,
            demuxed:    true,
            barcode:    barcode,
            parent:     meta.parent,
            parent_bam: meta.parent
        ], f)
    }.findAll { it != null }
}

// ---------------------------------------------------------------------------
// Workflow
// ---------------------------------------------------------------------------

workflow DEMUX {
    take:
        samplesheet_ch   // channel of rows (Map) parsed from samples.csv

    main:
        // ---- 1. Split rows into demux vs. pass-through paths ---------------
        rows = samplesheet_ch
            .map { row ->
                def needs_demux = (row.demux ?: 'false').toString().toLowerCase() in ['true', 't', '1', 'yes']
                row + [demux: needs_demux, reads: readsForRow(row)]
            }

        no_demux_rows = rows.filter { !it.demux }
        demux_rows    = rows.filter {  it.demux }

        // ---- 2. Pass-through path: emit (meta, reads) directly -------------
        passthrough_ch = no_demux_rows.map { row ->
            def meta = [
                id:         row.sample_id,
                sample_id:  row.sample_id,
                patient:    row.patient_id,
                timepoint:  row.timepoint,
                demuxed:    false,
                barcode:    null,
                parent:     row.reads.getName(),
                parent_bam: row.reads.getName()
            ]
            tuple(meta, row.reads)
        }

        passthrough_bam_ch   = passthrough_ch.filter { meta, f -> !isFastqFile(f) }
        passthrough_fastq_ch = passthrough_ch.filter { meta, f ->  isFastqFile(f) }

        // ---- 3. Demux path: group by parent reads file, run lima once ------
        // A "group" = all final samples that share the same parent reads file
        // and barcode FASTA. We collect the (barcode, sample_id) pairs
        // alongside so we can route per-barcode outputs back afterwards.
        grouped = demux_rows
            .map { row -> tuple(
                row.reads.toAbsolutePath().toString() + '|' +
                file(row.barcode_fasta).toAbsolutePath().toString(),
                row
            ) }
            .groupTuple()
            .map { key, rows_ ->
                def parent_path   = key.split('\\|')[0]
                def barcodes_path = key.split('\\|')[1]
                def parent_file   = file(parent_path)
                def parent_id     = parent_file.getName()
                                        .replaceAll(/(?i)\.(bam|fastq|fq|fasta|fa)(\.gz)?$/, '')

                // '*' in the barcode column means "take every barcode lima finds".
                def auto_rows = rows_.findAll { it.barcode?.toString()?.trim() == '*' }
                def auto_route = auto_rows.size() > 0
                if (auto_route && auto_rows.size() != rows_.size()) {
                    log.warn "[DEMUX] ${parent_file.getName()} mixes barcode='*' with explicit " +
                             "barcodes; taking ALL barcodes and ignoring the explicit rows."
                }

                def meta = [
                    id:           "demux_${parent_id}",
                    parent:       parent_file.getName(),
                    parent_bam:   parent_file.getName(),   // kept for back-compat
                    parent_id:    parent_id,
                    is_fastq:     isFastqFile(parent_file),
                    auto_route:   auto_route,
                    n_samples:    auto_route ? 0 : rows_.size(),
                    // Routing table: barcode -> sample_id (used after lima)
                    routing:      auto_route ? [:] : rows_.collectEntries { r ->
                                     [(r.barcode.toString()): r.sample_id]
                                  },
                    // Carry per-sample metadata for later annotation
                    samples_meta: auto_route ? [:] : rows_.collectEntries { r -> [(r.sample_id): [
                                     patient_id: r.patient_id,
                                     timepoint:  r.timepoint
                                  ]] }
                ]
                tuple(meta, parent_file, file(barcodes_path))
            }

        LIMA(grouped)

        // ---- 4. Re-emit one (meta, reads) tuple per FINAL sample -----------
        // BAM parents land in LIMA.out.bam; FASTQ parents land in
        // LIMA.out.fastqgz (gzipped input) or LIMA.out.fastq (plain).
        demuxed_bam_ch = LIMA.out.bam
            .flatMap { meta, files -> routeLimaOutputs(meta, files) }

        demuxed_fastq_ch = LIMA.out.fastqgz
            .mix(LIMA.out.fastq)
            .flatMap { meta, files -> routeLimaOutputs(meta, files) }

        // ---- 5. Concatenate both paths into per-sample channels ------------
        per_sample_bam   = passthrough_bam_ch.mix(demuxed_bam_ch)
        per_sample_fastq = passthrough_fastq_ch.mix(demuxed_fastq_ch)

    emit:
        bam     = per_sample_bam       // tuple(meta, bam)   for downstream stages
        fastq   = per_sample_fastq     // tuple(meta, fastq) for downstream stages
        summary = LIMA.out.summary     // optional QC artefacts
        counts  = LIMA.out.counts
        report  = LIMA.out.report
}
