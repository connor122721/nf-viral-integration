// modules/demux.nf
// ============================================================================
// Optional lima-based demultiplexing of multiplexed PacBio HiFi BAMs.
//
// Samplesheet schema (one row per FINAL sample after demultiplexing):
//
//   sample_id, patient_id, timepoint, bam, demux, barcode, barcode_fasta
//   PT1_T0,    PT1,        T0,        run1.bam, true,  bc1001--bc1001, refs/bc.fa
//   PT1_T1,    PT1,        T1,        run1.bam, true,  bc1002--bc1002, refs/bc.fa
//   PT2_T0,    PT2,        T0,        run2.bam, false, ,
//
//   - Rows with demux == 'true'  share a `bam` and run lima ONCE on it; each
//     row's `barcode` selects the corresponding lima output BAM.
//   - Rows with demux == 'false' (or missing/empty) bypass lima entirely.
//
// Public lima module is vendored at modules/nf-core/lima/main.nf.
// ============================================================================

include { LIMA } from './nf-core/lima/main'

workflow DEMUX {
    take:
        samplesheet_ch   // channel of rows (Map) parsed from samples.csv

    main:
        // ---- 1. Split rows into demux vs. pass-through paths ---------------
        rows = samplesheet_ch
            .map { row ->
                def needs_demux = (row.demux ?: 'false').toString().toLowerCase() in ['true', 't', '1', 'yes']
                row + [demux: needs_demux]
            }

        no_demux_rows = rows.filter { !it.demux }
        demux_rows    = rows.filter {  it.demux }

        // ---- 2. Pass-through path: emit (meta, bam) directly ---------------
        passthrough_ch = no_demux_rows.map { row ->
            def meta = [
                id:        row.sample_id,
                sample_id: row.sample_id,
                patient:   row.patient_id,
                timepoint: row.timepoint,
                demuxed:   false,
                parent_bam: file(row.bam).getName()
            ]
            tuple(meta, file(row.bam))
        }

        // ---- 3. Demux path: group by parent bam, run lima once per group ---
        // A "group" = all final samples that share the same parent BAM and
        // barcode FASTA. We collect the list of (barcode, sample_id) pairs
        // alongside so we can route per-barcode outputs back afterwards.
        grouped = demux_rows
            .map { row -> tuple(
                file(row.bam).toAbsolutePath().toString() + '|' +
                file(row.barcode_fasta).toAbsolutePath().toString(),
                row
            ) }
            .groupTuple()
            .map { key, rows_ ->
                def parent_bam_path = key.split('\\|')[0]
                def barcodes_path   = key.split('\\|')[1]
                def parent_id       = file(parent_bam_path).getBaseName()
                def meta = [
                    id:           "demux_${parent_id}",
                    parent_bam:   file(parent_bam_path).getName(),
                    n_samples:    rows_.size(),
                    // Routing table: barcode -> sample_id (used after lima)
                    routing:      rows_.collectEntries { r ->
                                     [(r.barcode.toString()): r.sample_id]
                                  },
                    // Carry per-sample metadata for later annotation
                    samples_meta: rows_.collectEntries { r -> [(r.sample_id): [
                                     patient_id: r.patient_id,
                                     timepoint:  r.timepoint
                                  ]] }
                ]
                tuple(meta, file(parent_bam_path), file(barcodes_path))
            }

        LIMA(grouped)

        // ---- 4. Re-emit one (meta, bam) tuple per FINAL sample -------------
        // lima writes <prefix>.<barcode>.bam (e.g. demux_run1.bc1001--bc1001.bam).
        // We pair each output BAM with the right sample_id via meta.routing.
        per_sample_ch = LIMA.out.bam
            .flatMap { meta, bams ->
                def bam_list = bams instanceof List ? bams : [bams]
                bam_list.collect { b ->
                    // Extract the barcode token from the filename:
                    //   demux_run1.bc1001--bc1001.bam -> bc1001--bc1001
                    def name    = b.getName()
                    def stem    = name.replaceFirst(/\.bam$/, '')
                    def barcode = stem.contains('.') ?
                                  stem.substring(stem.lastIndexOf('.') + 1) : stem
                    def sample_id = meta.routing[barcode]
                    if (!sample_id) {
                        log.warn "[DEMUX] lima produced BAM for unmapped barcode '${barcode}' " +
                                 "from parent ${meta.parent_bam} — skipping. " +
                                 "(routing keys: ${meta.routing.keySet()})"
                        return null
                    }
                    def smeta = meta.samples_meta[sample_id]
                    def out_meta = [
                        id:         sample_id,
                        sample_id:  sample_id,
                        patient:    smeta?.patient_id,
                        timepoint:  smeta?.timepoint,
                        demuxed:    true,
                        barcode:    barcode,
                        parent_bam: meta.parent_bam
                    ]
                    tuple(out_meta, b)
                }.findAll { it != null }
            }

        // ---- 5. Concatenate both paths into a single per-sample channel ---
        per_sample_bam = passthrough_ch.mix(per_sample_ch)

    emit:
        bam     = per_sample_bam       // tuple(meta, bam) for downstream stages
        summary = LIMA.out.summary     // optional QC artefacts
        counts  = LIMA.out.counts
        report  = LIMA.out.report
}
