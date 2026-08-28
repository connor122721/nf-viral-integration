#!/bin/env nextflow

/*
========================================================================================
    Viral Integration Detection Pipeline 
========================================================================================
    By: Connor S. Murray, PhD
    Based on: SMRTCap methodology (Smith Lab, University of Louisville)

    Workflow:
      0. [ Optional ] Lima demultiplexing of multiplexed HiFi BAMs.
      1. Select best viral reference via competitive mapping.
      2. Identify integration sites in human genome.
      3. Annotate integration sites and assign persistent clonal IDs and proviral integrity.
========================================================================================
*/
nextflow.enable.dsl = 2

// HELP MESSAGE
def helpMessage() {
    log.info"""
    ========================================================================================
    VIRAL INTEGRATION DETECTION PIPELINE
    ========================================================================================

    Usage:
        nextflow run main.nf --host_genome GRCh38.fa --viral_genomes "panel/*.fa" --annotation host.gtf

    Required Arguments:
        --host_genome           Path to host genome FASTA
        --viral_genomes         Glob pattern for viral reference panel (e.g., "hiv_refs/*.fa")

    Host Genome Annotation:
        --annotation            Host genome annotation file.
                                Accepted formats (plain or .gz): .gtf .gff .gff3

    Input Options (choose one):
        --samples               Path to samplesheet CSV.
                                Required column: sample
                                Reads:      'fastq' OR 'bam' (or both)
                                Tracking:   'patient_id', 'timepoint'  (optional but
                                            required if you want clonal-ID tracking
                                            across time-course samples)
                                Demux:      'demux' (true/false, default false),
                                            'barcode' (lima output token, e.g.
                                            'bc1001--bc1001'), 'barcode_fasta'
                                            (path to barcode FASTA)
                                            Required when demux == true.

                                Each row needs EITHER 'fastq' OR 'bam' populated.
                                Leave the unused cell empty. FASTQ rows always
                                bypass lima; BAM rows follow the demux settings.

                                Example mixing demuxed + multiplexed + FASTQ inputs:
                                  sample,patient_id,timepoint,fastq,bam,demux,barcode,barcode_fasta
                                  PT01_T0,PT01,T0,,runs/run1.bam,true,bc1001--bc1001,refs/bc.fa
                                  PT01_T1,PT01,T1,,runs/run1.bam,true,bc1002--bc1002,refs/bc.fa
                                  PT02_T0,PT02,T0,,runs/already.demuxed.bam,false,,
                                  STC1644,STC,T0,/data/STC1644.fastq.gz,,false,,

                                FASTQ-only sheets are fine — the 'bam' column can
                                be blank or absent entirely:
                                  sample,patient_id,timepoint,fastq
                                  STC1644,STC,T0,/data/STC1644.fastq.gz
                                  STC1649,STC,T1,/data/STC1649.fastq.gz

                                Multiplexed FASTQ demultiplexed by lima — set
                                demux=true and give the barcode token. Rows that
                                share a fastq + barcode_fasta run lima once:
                                  sample,patient_id,timepoint,fastq,demux,barcode,barcode_fasta
                                  PT03_T0,PT03,T0,runs/pool.fastq.gz,true,bc1003--bc1003,refs/bc.fa
                                  PT03_T1,PT03,T1,runs/pool.fastq.gz,true,bc1004--bc1004,refs/bc.fa

                                Take every barcode from a pool without naming
                                them up front:
                                  sample,fastq,demux,barcode,barcode_fasta
                                  pool1,runs/pool.fastq.gz,true,*,refs/bc.fa

        --patient_dir           Directory containing patient BAM/FASTQ files
        --patient_bam           Single patient BAM file
        --patient_fastq         Single patient FASTQ file (already demultiplexed)
        --multiplexed_fastq     Multiplexed FASTQ(s) to demultiplex with lima.
                                Glob allowed. Requires --default_barcode_fasta.
                                Every barcode lima finds becomes a sample, named
                                <run>.<barcode> (e.g. run1.bc1001--bc1001).

    Demultiplexing (lima — modules/nf-core/lima/main.nf):
        --skip_demux            Bypass lima entirely; treat every row as already
                                demultiplexed regardless of its 'demux' value.
        --lima_args             String passed verbatim to lima.
                                Default: "--hifi-preset SYMMETRIC --min-score 80"
                                Note: --split-named is appended automatically
                                unless you supply your own --split* flag.
        --default_barcode_fasta Fallback barcode FASTA if a samplesheet row sets
                                demux=true but omits 'barcode_fasta'. Also the
                                barcode FASTA used by --multiplexed_fastq.

        FASTQ vs BAM demultiplexing: both are handled by the same lima module
        (modules/demux.nf). lima mirrors the input extension, so a BAM pool
        yields <prefix>.<barcode>.bam and a FASTQ pool yields
        <prefix>.<barcode>.fastq[.gz]. lima cannot emit BAM from FASTQ input,
        so demuxed FASTQs never enter BAM_TO_FASTQ.

        IMPORTANT: lima only writes one file per barcode when --split-named is
        in effect. That flag comes from the LIMA module's ext.args in
        nextflow.config, not from this pipeline, so make sure it is set.

        Set 'barcode' to '*' on a demux=true row to take EVERY barcode lima
        finds in that pool; samples are then named <pool>.<barcode>.

    Post-processing:
        --skip_repeatmasker     Skip RepeatMasker download + overlap [default: false]

    Output:
        --outdir                Output directory [default: ./output]

    Profiles:
        -profile slurm,singularity    Use Slurm scheduler with Apptainer/Singularity

    Examples:
        # Multi-sample run with mixed demuxed + multiplexed inputs
        nextflow run main.nf \\
            --samples samples.csv \\
            --host_genome t2t.fa \\
            --annotation chm13v2.0.gtf \\
            --viral_genomes "hiv_refs/*.fa" \\
            -profile slurm,singularity

        # Same, but skip lima (everything is pre-demultiplexed)
        nextflow run main.nf --samples samples.csv --skip_demux ...

    ========================================================================================
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

// ========================================================================================
// INPUT VALIDATION
// ========================================================================================

// True only if the samplesheet row has that column AND it is non-null / non-blank.
// Handles missing columns, empty cells ("PT02,,,"), and whitespace-only cells.
def hasVal(row, col) {
    row instanceof Map && row.containsKey(col) && row[col] != null && row[col].toString().trim() != ''
}

// Clean up samplesheet header keys and values before anything downstream sees them.
//   - strips a UTF-8 BOM from the first column name  ("\uFEFFsample" -> "sample")
//   - strips stray whitespace / carriage returns from keys and values
//   - lowercases keys so 'Sample' / 'BAM' / 'FastQ' all work
//   - accepts common aliases for the sample-name column
def normaliseRow(row) {
    def clean = [:]
    row.each { k, v ->
        def key = k?.toString()?.replaceAll(/\uFEFF/, '')?.replaceAll(/\r/, '')?.trim()?.toLowerCase()
        if (key) clean[key] = (v instanceof CharSequence) ? v.toString().replaceAll(/\r/, '').trim() : v
    }

    if (!clean.sample) {
        def alias = ['sample_id','sampleid','sample_name','samplename','id','name']
                        .find { clean[it] != null && clean[it].toString().trim() != '' }
        if (alias) clean.sample = clean[alias]
    }

    if (!clean.sample) {
        throw new IllegalArgumentException(
            "Samplesheet row has no usable sample name.\n" +
            "       Columns found: ${clean.keySet().join(', ')}\n" +
            "       Expected a 'sample' column (aliases: sample_id, sample_name, id, name).\n" +
            "       If the header looks correct, the file likely has a UTF-8 BOM or CRLF\n" +
            "       line endings — check with:  head -1 samples.csv | xxd | head -2")
    }

    // modules/demux.nf reads row.sample_id; the rest of main.nf reads
    // row.sample. Keep both populated and identical so neither goes null.
    clean.sample_id = clean.sample
    clean
}

// Loose boolean parser for samplesheet cells ("true"/"T"/"1"/"yes"/"Y").
def isTrue(v) {
    v != null && v.toString().trim().toLowerCase() in ['true','t','1','yes','y']
}

// Fill in the fallback barcode FASTA and check that demux rows are complete.
// Applied to every row; rows that do not request demultiplexing pass untouched.
def resolveDemuxRow(row) {
    if (!isTrue(row.demux)) return row

    if (!hasVal(row, 'barcode_fasta') && params.default_barcode_fasta) {
        row.barcode_fasta = params.default_barcode_fasta
    }
    if (!hasVal(row, 'barcode_fasta')) {
        throw new IllegalArgumentException(
            "Samplesheet row '${row.sample}' sets demux=true but has no 'barcode_fasta'.\n" +
            "       Add the column or pass --default_barcode_fasta.")
    }
    if (!hasVal(row, 'barcode')) {
        throw new IllegalArgumentException(
            "Samplesheet row '${row.sample}' sets demux=true but has no 'barcode'.\n" +
            "       Expected a lima --split-named token (e.g. 'bc1001--bc1001'),\n" +
            "       or '*' to take every barcode lima finds in that pool.")
    }
    row
}

def validateInputs() {
    def errors = []

    if (!params.host_genome) {
        errors << "ERROR: --host_genome is required."
    } else {
        def genome = params.host_genome.toString()
        if (!(genome ==~ /.*\.(fa|fasta|fna)(\.gz)?$/)) {
            errors << "ERROR: --host_genome '${genome}' has an unrecognised extension.\n" +
                      "       Accepted: .fa  .fasta  .fa.gz  .fasta.gz  .fna  .fna.gz"
        }
    }

    if (params.annotation) {
        def annot = params.annotation.toString()
        if (!(annot ==~ /.*\.(gtf|gff|gff3)(\.gz)?$/)) {
            errors << "ERROR: --annotation '${annot}' has an unrecognised extension.\n" +
                      "       Accepted: .gtf  .gff  .gff3  .gtf.gz  .gff.gz  .gff3.gz"
        }
    }

    if (params.report_genome && !(params.report_genome.toString().toLowerCase() in ['hg38', 't2t'])) {
        errors << "ERROR: --report_genome must be 'hg38' or 't2t' (got '${params.report_genome}')."
    }

    if (params.html_mode && !(params.html_mode.toString().toLowerCase() in ['single','multi','both'])) {
        errors << "ERROR: --html_mode must be 'single', 'multi', or 'both' (got '${params.html_mode}')."
    }

    if (errors) {
        log.error "\n" + errors.join("\n")
        exit 1
    }
}

validateInputs()

// ========================================================================================
// PRINT PARAMETER SUMMARY
// ========================================================================================

log.info ""
log.info "========================================================================================"
log.info "HIV SMRTCap INTEGRATION DETECTION PIPELINE"
log.info "========================================================================================"
log.info "Host genome       : ${params.host_genome}"
log.info "Viral genomes     : ${params.viral_genomes ?: params.viral_genome}"
log.info "Host annotation   : ${params.annotation}"
log.info "Demultiplex (lima): ${params.skip_demux ? 'SKIPPED' : 'enabled (per-row demux flag; BAM + FASTQ)'}"
log.info "Barcode FASTA     : ${params.default_barcode_fasta ?: 'per-row only'}"
log.info "Report genome     : ${params.report_genome}"
log.info "Output directory  : ${params.outdir}"
log.info "========================================================================================"
log.info ""

// ========================================================================================
//  PROCESSES
// ========================================================================================

// Convert BAM to FASTQ
process BAM_TO_FASTQ {
    tag "${bam_file.baseName}"
    publishDir "${params.outdir}/converted_fastq", mode: 'link'
    container params.container

    input:
        path bam_file

    output:
        path "*.fastq.gz", emit: fastq

    script:
    """
    samtools fastq -@ ${params.threads} ${bam_file} > ${bam_file.baseName}.fastq
    pigz ${bam_file.baseName}.fastq
    """
}

// Map reads to each viral reference genome (competitive mapping)
process MULTI_REFERENCE_MAPPING {
    tag "${sample_id}_vs_${viral_genome.baseName}"
    publishDir "${params.outdir}/01_reference_selection/${sample_id}", mode: 'link'
    container params.container

    input:
        tuple val(sample_id), path(reads)
        each path(viral_genome)
        path(host_genome)

    output:
        tuple val(sample_id), path(viral_genome), path("*sorted.sam"), path("*primary.sam"), path("*stats.txt"),
             path("*sorted.sam.fa"), path("*viralreads.fa"), path("*hostflanks.fa"), emit: results_selectBest
        tuple val(sample_id), path("*.pbmarkdup.log")
        tuple val(sample_id), path("*.readnames.txt")
        tuple val(sample_id), path("*.sorted.bam"), emit: sorted_bam
        path("*sorted.sam.fa"), emit: mask_fasta
        path("*viralreads.fa"), emit: unmask_fasta
        path("*viralhits*.fastq.gz")
        path("*viralhits.sam")

    script:
        def ref_name = viral_genome.baseName
        def sample_id_i = sample_id.toString().replaceAll(/\.bam$/, '')
    """
    # Build combined host+viral reference once
    mkdir -p ${projectDir}/tmp
    if [ ! -f ${projectDir}/tmp/${ref_name}_hybridhost.fa ]; then
        cat ${viral_genome} ${host_genome} > ${projectDir}/tmp/${ref_name}_hybridhost.fa
    fi

    # Pass 1: select reads that hit the virus genome
    minimap2 -t ${params.threads} -m 0 -Y -ax map-hifi --score-N=0 --sam-hit-only ${viral_genome} ${reads} > \\
        ${sample_id_i}_vs_${ref_name}_viralhits.sam
    
    samtools fastq -@ ${params.threads} -F 0x904 ${sample_id_i}_vs_${ref_name}_viralhits.sam | gzip > \\
        ${sample_id_i}_viralhits.fastq.gz

    # Get viral reads from alignment
    zcat ${sample_id_i}_viralhits.fastq.gz | grep "@" | sed 's/@//g' | cut -f1 -d" " > \\
        ${sample_id_i}_vs_${ref_name}.viralhits.readnames.txt

    # Mask viral reads
    python ${projectDir}/bin/mask.py ${sample_id_i}_vs_${ref_name}_viralhits.sam > \\
        ${sample_id_i}_vs_${ref_name}.sorted.sam.fa

    # Extract non-viral flanks - this is host and consolidates any reads with N gaps
    python ${projectDir}/bin/get_flanks.py ${sample_id_i}_vs_${ref_name}.sorted.sam.fa > \\
        ${sample_id_i}_vs_${ref_name}.hostflanks.fa

    # Extract viral reads only
    python ${projectDir}/bin/unmask.py \\
        ${sample_id_i}_vs_${ref_name}_viralhits.sam \\
        ${sample_id_i}_vs_${ref_name}.sorted.sam.fa > \\
        ${sample_id_i}_vs_${ref_name}.viralreads.fa

    # Pass 2: map selected reads to the combined reference
    minimap2 -t ${params.threads} -m 0 -Y -ax map-hifi --score-N=0 \\
        ${projectDir}/tmp/${ref_name}_hybridhost.fa ${sample_id_i}_viralhits.fastq.gz | \\
        samtools view -h -F 4 -b | \\
        samtools sort -@ ${params.threads} -o ${sample_id_i}_vs_${ref_name}.sorted.bam

    # Index required by parseSAM's region query (samtools view bam HIV:1-1000000)
    samtools index ${sample_id_i}_vs_${ref_name}.sorted.bam
    samtools view ${sample_id_i}_vs_${ref_name}.sorted.bam | cut -f1 | sort -u > \\
        ${sample_id_i}_vs_${ref_name}.allmapped.readnames.txt

    # pbmarkdup input
    samtools fastq -@ ${params.threads} ${sample_id_i}_vs_${ref_name}.sorted.bam > \\
        ${sample_id_i}_vs_${ref_name}.fastq

    pbmarkdup ${sample_id_i}_vs_${ref_name}.fastq \\
        ${sample_id_i}_vs_${ref_name}.markdup.fastq \\
        --dup-file ${sample_id_i}_vs_${ref_name}.dups.fastq \\
        --log-level INFO > ${sample_id_i}_vs_${ref_name}.pbmarkdup.log

    grep "@" ${sample_id_i}_vs_${ref_name}.dups.fastq | sed 's/@//g' | cut -f1 -d" " > \\
        ${sample_id_i}_vs_${ref_name}.dups.readnames.txt

    # Full-flag SAM (supplementaries intact) for flagstat / inspection.
    # --qname-file ^${sample_id_i}_vs_${ref_name}.dups.readnames.txt
    samtools view -h ${sample_id_i}_vs_${ref_name}.sorted.bam \\
        -o ${sample_id_i}_vs_${ref_name}.sorted.sam

    samtools flagstat ${sample_id_i}_vs_${ref_name}.sorted.sam > \\
        ${sample_id_i}_vs_${ref_name}.stats.txt

    # Primary-only SAM for mask.py -> exactly one FASTA record per read.
    # --qname-file ^${sample_id_i}_vs_${ref_name}.dups.readnames.txt
    samtools view -h -F 0x900 ${sample_id_i}_vs_${ref_name}.sorted.bam \\
        -o ${sample_id_i}_vs_${ref_name}.primary.sam

     mv ${sample_id_i}_viralhits.fastq.gz ${sample_id_i}_${ref_name}_viralhits.fastq.gz
    """
}

// Select best viral reference based on alignment metrics
process SELECT_BEST_REFERENCE {
    tag "${sample_id}"
    publishDir "${params.outdir}/01_reference_selection/${sample_id}", mode: 'link'
    container params.container

    input:
        tuple val(sample_id), path(viral_genomes), path(sam_files_supp), 
            path(sam_files_nonsupp), path(stats_files), path(sort_fa), path(unmask_fa), path(host_flanks_fa)

    output:
        tuple val(sample_id), path("*best_reference.txt"), emit: best_ref_name
        tuple val(sample_id), path("*best_reference.fa"), emit: best_ref_fa
        tuple val(sample_id), path("*best_reference.sam"), emit: best_sam
        path "*mapping_comparison.txt", emit: comparison
        path "*detailed_metrics.txt", emit: metrics
        tuple val(sample_id), path("*final.masked.fa"), emit: mask_fa
        tuple val(sample_id), path("*final.unmasked.fa"), emit: unmask_fa
        tuple val(sample_id), path("*final.hostflanks.fa"), emit: host_flanks_fa

    script:
        def sam_list = sam_files_supp.collect{ it }.join(' ')
        def stats_list = stats_files.collect{ it }.join(' ')
        def genome_list = viral_genomes.collect{ it }.join(' ')
        def sample_id_i = sample_id.toString().replaceAll(/\.bam$/, '')
    """
    python ${projectDir}/bin/select_best_reference.py \\
        --sam-files ${sam_list} \\
        --stats-files ${stats_list} \\
        --viral-genomes ${genome_list} \\
        --output-best ${sample_id_i}_best_reference.txt \\
        --output-fa ${sample_id_i}_best_reference.fa \\
        --output-comparison ${sample_id_i}_mapping_comparison.txt \\
        --output-detailed ${sample_id_i}_detailed_metrics.txt

    BEST_REF=\$(head -n1 ${sample_id_i}_best_reference.txt )
    BEST_REF_vir=\$(head -n1 ${sample_id_i}_best_reference.txt | sed 's/.sorted//g')
    cp ${sample_id_i}_vs_\${BEST_REF}.sam ${sample_id_i}_vs_\${BEST_REF}_best_reference.sam
    cp ${sample_id_i}_vs_\${BEST_REF}.sam.fa ${sample_id_i}.final.masked.fa
    cp ${sample_id_i}_vs_\${BEST_REF_vir}.viralreads.fa ${sample_id_i}.final.unmasked.fa
    cp ${sample_id_i}_vs_\${BEST_REF_vir}.hostflanks.fa ${sample_id_i}.final.hostflanks.fa
    """
}

// Rename a demuxed FASTQ to <sample_id>.fastq.gz for QC + downstream
process RENAME_DEMUXED_FASTQ {
    tag "${meta.id}"
    container params.container

    input:
        tuple val(meta), path(fastq, stageAs: 'input/*')

    output:
        tuple val(meta), path("${meta.id}.fastq.gz"), emit: fastq

    script:
    """
    case "${fastq}" in
        *.gz)
            ln -sf "\$(readlink -f ${fastq})" ${meta.id}.fastq.gz
            ;;
        *)
            pigz -p ${params.threads} -c ${fastq} > ${meta.id}.fastq.gz
            ;;
    esac
    """
}

// Rename a demuxed BAM to <sample_id>.bam so downstream BAM_TO_FASTQ + sample_id
process RENAME_DEMUXED_BAM {
    tag "${meta.id}"
    container params.container

    input:
        tuple val(meta), path(bam, stageAs: 'input/*')

    output:
        path "${meta.id}.bam", emit: bam

    script:
    """
    ln -sf "\$(readlink -f ${bam})" ${meta.id}.bam
    """
}

// ========================================================================================
//  MODULE INCLUDES
// ========================================================================================

// QC modules
include { FASTQC } from './bin/qc_mods.nf'
include { QUALIMAP } from './bin/qc_mods.nf'
include { MULTIQC } from './bin/qc_mods.nf'

// Reference + downstream genomic processes
include { DEMUX } from './modules/demux'
include { POST_PROCESS } from './modules/post_process'
include { GFFCONVERT } from './bin/genomic_processes.nf'
include { UNMASK_SEQUENCES } from './bin/genomic_processes.nf'
include { EXTRACT_FLANKS } from './bin/genomic_processes.nf'
include { MAP_FLANKS_TO_HOST } from './bin/genomic_processes.nf'
include { CONFIRM_HOST_ALIGNMENTS } from './bin/genomic_processes.nf'
include { COMBINE_RESULTS } from './bin/genomic_processes.nf'
include { INTEGRATION_ANNOTATE } from './bin/integration_annotation.nf'
include { CREATE_HTML_REPORT } from './bin/integration_annotation.nf'

// ========================================================================================
// WORKFLOW
// ========================================================================================
workflow {

    // ==================================================================================
    // SETUP
    // ==================================================================================
    host_genome_ch = Channel.fromPath(params.host_genome, checkIfExists: true)
    viral_genomes_ch = Channel.fromPath(params.viral_genomes, checkIfExists: true)
    annotation_ch = params.annotation ? Channel.fromPath(params.annotation, checkIfExists: true) : Channel.empty()
    n_viral_refs = files(params.viral_genomes).size()

    // Script paths
    script_dir = "${projectDir}/bin"
    unmask_script_ch = Channel.fromPath("${script_dir}/unmask.py", checkIfExists: true)
    get_flanks_script_ch = Channel.fromPath("${script_dir}/get_flanks.py", checkIfExists: true)
    annotate_script_ch = Channel.fromPath("${script_dir}/1.annotate_integrations_HIV.R", checkIfExists:true)
    blast_script_ch = Channel.fromPath("${script_dir}/findHIVSIVGeneRegionsV9.6", checkIfExists: true)
    clone_calling_script_ch = Channel.fromPath("${script_dir}/parseSAM_v2.pl", checkIfExists: true)
    report_script_ch = Channel.fromPath("${script_dir}/3.Create_Sample_HTML.R", checkIfExists: true)

    // ==================================================================================
    // STEP 0: Prepare input reads
    //         - samplesheet / --multiplexed_fastq: optional lima demultiplex
    //           (BAM *and* FASTQ), BAM->FASTQ conversion, FASTQ pass-through
    //         - patient_dir / patient_bam / patient_fastq
    // ==================================================================================
    if (params.samples || params.multiplexed_fastq || params.patient_dir || params.patient_bam || params.patient_fastq) {

        // ROW-BASED MODES (samplesheet or --multiplexed_fastq) ------------------
        if (params.samples || params.multiplexed_fastq) {

            if (params.samples) {
                rows_ch = Channel
                    .fromPath(params.samples, checkIfExists: true)
                    .splitCsv(header: true, strip: true)
                    .map { raw ->
                        // Normalise header keys/values (BOM, CRLF, case, aliases) first.
                        def row = normaliseRow(raw)

                        // Every row must carry reads in at least one of the two columns.
                        if (!hasVal(row, 'fastq') && !hasVal(row, 'bam')) {
                            throw new IllegalArgumentException(
                                "Samplesheet row '${row.sample}' has neither a 'fastq' nor a " +
                                "'bam' value. Populate one of them.")
                        }
                        row
                    }

            } else {
                // --multiplexed_fastq: synthesise one demux row per input pool.
                // barcode='*' tells DEMUX to emit every barcode lima finds,
                // naming samples <pool>.<barcode>.
                if (!params.default_barcode_fasta) {
                    throw new IllegalArgumentException(
                        "--multiplexed_fastq requires --default_barcode_fasta (the barcode FASTA).")
                }
                if (params.skip_demux) {
                    throw new IllegalArgumentException(
                        "--multiplexed_fastq is inherently multiplexed; it cannot be combined " +
                        "with --skip_demux. Use --patient_fastq for pre-demultiplexed reads.")
                }

                rows_ch = Channel
                    .fromPath(params.multiplexed_fastq, checkIfExists: true)
                    .map { f ->
                        def run_id = f.name.replaceAll(/\.(fastq|fq)(\.gz)?$/, '')
                        [ sample:        run_id,
                          sample_id:     run_id,
                          patient_id:    null,
                          timepoint:     null,
                          fastq:         f.toString(),
                          demux:         'true',
                          barcode:       '*',
                          barcode_fasta: params.default_barcode_fasta ]
                    }
            }

            if (!params.skip_demux) {

                // DEMUX handles all four cases in one pass: demux/pass-through
                // x BAM/FASTQ. It emits tuple(meta, reads) on .bam and .fastq.
                DEMUX(rows_ch.map { row -> resolveDemuxRow(row) })

                // --- BAM side: rename to <sample>.bam, then convert to FASTQ ---
                RENAME_DEMUXED_BAM(DEMUX.out.bam)
                BAM_TO_FASTQ(RENAME_DEMUXED_BAM.out.bam)
                bam_converted_ch = BAM_TO_FASTQ.out.fastq
                    .map { f -> tuple(f.baseName.replaceAll(/\.(fastq|fq)(\.gz)?$/, ''), f) }

                // --- FASTQ side: rename/compress to <sample>.fastq.gz ----------
                // lima cannot emit BAM from FASTQ input, so these never enter
                // BAM_TO_FASTQ. If the sheet is FASTQ-only, DEMUX.out.bam is
                // empty and the BAM branch above simply spawns zero tasks.
                RENAME_DEMUXED_FASTQ(DEMUX.out.fastq)
                fastq_ready_ch = RENAME_DEMUXED_FASTQ.out.fastq
                    .map { meta, f -> tuple(meta.id, f) }

                input_reads_ch = fastq_ready_ch.mix(bam_converted_ch)

            } else {
                // --skip_demux: bypass lima entirely, mirroring DEMUX's
                // precedence rule that 'bam' wins when a row sets both.
                bam_files_ch = rows_ch
                    .filter { row -> hasVal(row, 'bam') }
                    .map    { row -> file(row.bam, checkIfExists: true) }

                BAM_TO_FASTQ(bam_files_ch)
                bam_converted_ch = BAM_TO_FASTQ.out.fastq
                    .map { f -> tuple(f.baseName.replaceAll(/\.(fastq|fq)(\.gz)?$/, ''), f) }

                fastq_ready_ch = rows_ch
                    .filter { row -> hasVal(row, 'fastq') && !hasVal(row, 'bam') }
                    .map    { row -> tuple(row.sample, file(row.fastq, checkIfExists: true)) }

                input_reads_ch = fastq_ready_ch.mix(bam_converted_ch)
            }


        } else if (params.patient_dir) {

            bam_files_ch = Channel.fromPath("${params.patient_dir}/**/*.bam", checkIfExists: false)
                .mix(Channel.fromPath("${params.patient_dir}/*.bam", checkIfExists: false))
                .filter { it.exists() }
                .filter { !it.name.toLowerCase().contains('unassigned') }
                .filter { !it.toUriString().toLowerCase().contains('fail_reads') }
                .unique()

            fastq_files_ch = Channel.fromPath("${params.patient_dir}/**/*.{fastq,fq,fastq.gz,fq.gz}", checkIfExists: false)
                .mix(Channel.fromPath("${params.patient_dir}/*.{fastq,fq,fastq.gz,fq.gz}", checkIfExists: false))
                .filter { it.exists() }
                .filter { !it.name.toLowerCase().contains('unassigned') }
                .unique()

            BAM_TO_FASTQ(bam_files_ch)
            input_reads_ch = BAM_TO_FASTQ.out.fastq.mix(fastq_files_ch)
                .map { file -> def sample_id = file.baseName.replaceAll(/\.(fastq|fq)(\.gz)?$/, '')
                    tuple(sample_id, file) }

        } else if (params.patient_bam) {
            bam_ch = Channel.fromPath(params.patient_bam, checkIfExists: true)
            BAM_TO_FASTQ(bam_ch)
            input_reads_ch = BAM_TO_FASTQ.out.fastq
                .map { file -> def sample_id = file.baseName.replaceAll(/\.(fastq|fq)(\.gz)?$/, '')
                    tuple(sample_id, file) }
        } else {
            input_reads_ch = Channel.fromPath(params.patient_fastq, checkIfExists: true)
                .map { file -> def sample_id = file.baseName.replaceAll(/\.(fastq|fq)(\.gz)?$/, '')
                    tuple(sample_id, file) }
        }
    }

    // Convert/pass-through annotation file to GTF regardless of input format
    GFFCONVERT(annotation_ch)
    gtf_ch = GFFCONVERT.out.gtf.first()
    FASTQC(input_reads_ch)

    // ==================================================================================
    // STEP 1: Select best viral reference via competitive mapping/chimeric per sample
    // ==================================================================================
    MULTI_REFERENCE_MAPPING(input_reads_ch, viral_genomes_ch, host_genome_ch.first())

    grouped_results = MULTI_REFERENCE_MAPPING.out.results_selectBest
        .groupTuple(by: 0, size: n_viral_refs, remainder: true)
        .map { sample_id, viral_genomes, sams, primary_sams, stats, fasta, unmask_fa, host_flank ->
            tuple(sample_id, viral_genomes, sams, primary_sams,stats, fasta, unmask_fa, host_flank) }

    //QUALIMAP(MULTI_REFERENCE_MAPPING.out.sorted_bam)
    SELECT_BEST_REFERENCE(grouped_results)

    // every sample that entered the step
    all_ids = grouped_results.map { tuple(it[0], 'seen') }
    ok_ids = SELECT_BEST_REFERENCE.out.best_ref_name.map { tuple(it[0], 'ok') }

    failed_ids = all_ids
        .join(ok_ids, remainder: true)
        .filter { id, seen, ok -> ok == null }
        .map { id, seen, ok -> id }

    failed_ids
        .collectFile(name: 'failed_best_reference_samples.txt', newLine: true, sort: true,
                    storeDir: "${params.outdir}")

    annotate_input = SELECT_BEST_REFERENCE.out.best_ref_fa
        .join(SELECT_BEST_REFERENCE.out.unmask_fa)
        .join(SELECT_BEST_REFERENCE.out.host_flanks_fa)
        .join(SELECT_BEST_REFERENCE.out.best_sam)

    // Integrate the annotations
    INTEGRATION_ANNOTATE(annotate_input,
                         clone_calling_script_ch.first(),
                         annotate_script_ch.first(),
                         blast_script_ch.first(),
                         report_script_ch.first(),
                         gtf_ch,
			 params.repeatmasker_cache_dir)

    //CREATE_HTML_REPORT(INTEGRATION_ANNOTATE.out.csv_ann.collect(),
    //                   project_report_script_ch.first())

    // QC output
    MULTIQC(FASTQC.out.zip
        .mix(MULTI_REFERENCE_MAPPING.out.results_selectBest.map { it[3] })
        //.mix(QUALIMAP.out.results)
        .collect())
}

// ========================================================================================
// COMPLETION
// ========================================================================================
workflow.onComplete {
    log.info ""
    log.info "========================================================================================"
    log.info "VIRAL INTEGRATION PIPELINE COMPLETE"
    log.info "========================================================================================"
    log.info "Status:           ${workflow.success ? 'SUCCESS' : 'FAILED'}"
    log.info "Duration:         ${workflow.duration}"
    log.info "Output directory: ${params.outdir}"
    log.info ""
    log.info "Key Results:"
    log.info "  01_reference_selection/   - Best viral reference & initial mapping"
    log.info "  final_results/            - Integration sites in human genome"
    log.info ""
    log.info "========================================================================================"
}
