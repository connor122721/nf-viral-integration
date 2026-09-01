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

                                Example mixing demuxed + multiplexed inputs:
                                  sample,patient_id,timepoint,bam,demux,barcode,barcode_fasta
                                  PT01_T0,PT01,T0,runs/run1.bam,true,bc1001--bc1001,refs/bc.fa
                                  PT01_T1,PT01,T1,runs/run1.bam,true,bc1002--bc1002,refs/bc.fa
                                  PT02_T0,PT02,T0,runs/already.demuxed.bam,false,,
                                  STC1644,STC,T0,/data/STC1644.fastq.gz,false,,

        --patient_dir           Directory containing patient BAM/FASTQ files
        --patient_bam           Single patient BAM file
        --patient_fastq         Single patient FASTQ file

    Demultiplexing (lima — modules/nf-core/lima/main.nf):
        --skip_demux            Bypass lima entirely; treat every row as already
                                demultiplexed regardless of its 'demux' value.
        --lima_args             String passed verbatim to lima.
                                Default: "--hifi-preset SYMMETRIC --min-score 80"
        --default_barcode_fasta Fallback barcode FASTA if a samplesheet row sets
                                demux=true but omits 'barcode_fasta'.

    Post-processing (clonal IDs, repeats, circos, HTML reports):
        --skip_repeatmasker     Skip RepeatMasker download + overlap [default: false]
        --report_genome         Used by circos + repeats: hg38 or t2t [default: hg38]
        --html_mode             single | multi | both [default: both]

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
log.info "Demultiplex (lima): ${params.skip_demux ? 'SKIPPED' : 'enabled (per-row demux flag)'}"
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
        tuple val(sample_id), path(viral_genome), path("*sorted.sam"), path("*primary.sam"), path("*stats.txt"), path("*sorted.sam.fa"), emit: results_selectBest
        tuple val(sample_id), path("*.pbmarkdup.log")
        tuple val(sample_id), path("*.readnames.txt")
        tuple val(sample_id), path("*.sorted.bam"), emit: sorted_bam
        path("*sorted.sam.fa"), emit: unmask_fasta

    script:
        def ref_name = viral_genome.baseName
        def sample_id_i = sample_id.replaceAll(/.bam$/, '')
    """
    # Build combined host+viral reference once
    mkdir -p ${projectDir}/tmp
    if [ ! -f ${projectDir}/tmp/${ref_name}_hybridhost.fa ]; then
        cat ${viral_genome} ${host_genome} > ${projectDir}/tmp/${ref_name}_hybridhost.fa
    fi

    # Pass 1: select reads that hit virus
    minimap2 -t ${params.threads} -m 0 -Y -ax map-hifi --score-N=0 \\
        ${viral_genome} ${reads} | \\
        samtools fastq -@ ${params.threads} -F 0x904 - | \\
        gzip > ${sample_id_i}_viralhits.fastq.gz

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
        --log-file ${sample_id_i}_vs_${ref_name}.pbmarkdup.log

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

    python ${projectDir}/bin/mask.py ${sample_id_i}_vs_${ref_name}.primary.sam > \\
        ${sample_id_i}_vs_${ref_name}.sorted.sam.fa
    """
}

// Select best viral reference based on alignment metrics
process SELECT_BEST_REFERENCE {
    tag "${sample_id}"
    publishDir "${params.outdir}/01_reference_selection/${sample_id}", mode: 'link'
    container params.container

    input:
        tuple val(sample_id), path(viral_genomes), path(sam_files_supp), path(sam_files_nonsupp), path(stats_files), path(sort_fa)

    output:
        tuple val(sample_id), path("*best_reference.txt"), emit: best_ref_name
        tuple val(sample_id), path("*best_reference.fa"), emit: best_ref_fa
        tuple val(sample_id), path("*best_reference.sam"), emit: best_sam
        path "*mapping_comparison.txt", emit: comparison
        path "*detailed_metrics.txt", emit: metrics
        tuple val(sample_id), path("*final.masked.fa"), emit: mask_fa

    script:
        def sam_list = sam_files_supp.collect{ it }.join(' ')
        def stats_list = stats_files.collect{ it }.join(' ')
        def genome_list = viral_genomes.collect{ it }.join(' ')
        def sample_id_i = sample_id.replaceAll(/.bam$/, '')
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
    cp ${sample_id_i}_vs_\${BEST_REF}.sam ${sample_id_i}_vs_\${BEST_REF}_best_reference.sam
    cp ${sample_id_i}_vs_\${BEST_REF}.sam.fa ${sample_id_i}.final.masked.fa
    """
}

// Rename a demuxed BAM to <sample_id>.bam so downstream BAM_TO_FASTQ + sample_id
process RENAME_DEMUXED_BAM {
    tag "${meta.id}"
    container params.container

    input:
        tuple val(meta), path(bam)

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
    //         - samplesheet: optional lima demultiplex, BAM->FASTQ, FASTQ pass-through
    //         - patient_dir / patient_bam / patient_fastq
    // ==================================================================================
    if (params.samples || params.patient_dir || params.patient_bam || params.patient_fastq) {

        // SAMPLESHEET MODE -----------------------------------------------------
        if (params.samples) {

            samplesheet_ch = Channel
                .fromPath(params.samples, checkIfExists: true)
                .splitCsv(header: true, strip: true)

            // ---- Optional lima demultiplex on rows with demux=true ----------
            if (!params.skip_demux) {
                // Tag rows for the DEMUX sub-workflow. DEMUX returns
                // tuple(meta, demuxed_bam) for demux=true rows and the same
                // shape (meta, original_bam) for the pass-through path.
                demux_input_rows = samplesheet_ch.map { row ->
                    // Allow either a per-row barcode_fasta OR fall back to params.default_barcode_fasta
                    if (row.demux?.toString()?.toLowerCase() in ['true','t','1','yes']
                        && (!row.barcode_fasta || row.barcode_fasta.trim() == '')
                        && params.default_barcode_fasta) {
                        row.barcode_fasta = params.default_barcode_fasta
                    }
                    row
                }
                DEMUX(demux_input_rows)

                // Rename so each demuxed BAM is named ${sample}.bam (lima
                // emits <prefix>.<barcode>.bam by default).
                RENAME_DEMUXED_BAM(DEMUX.out.bam)

                // FASTQ rows from samplesheet pass straight through; demux
                // outputs and non-demux BAM rows enter BAM_TO_FASTQ.
                fastq_rows_ch = samplesheet_ch
                    .filter { row -> row.containsKey('fastq') && row.fastq && row.fastq.trim() != '' }
                    .map { row -> tuple(row.sample, file(row.fastq, checkIfExists: true)) }

                bam_files_ch = RENAME_DEMUXED_BAM.out.bam

                BAM_TO_FASTQ(bam_files_ch)
                bam_converted_ch = BAM_TO_FASTQ.out.fastq
                    .map { f -> tuple(f.baseName.replaceAll(/\.(fastq|fq)(\.gz)?$/, ''), f) }

                input_reads_ch = fastq_rows_ch.mix(bam_converted_ch)

            } else {
                // skip_demux
                bam_rows_ch = samplesheet_ch
                    .filter { row -> row.containsKey('bam') && row.bam && row.bam.trim() != '' }
                    .map    { row -> file(row.bam, checkIfExists: true) }

                BAM_TO_FASTQ(bam_rows_ch)
                bam_converted_ch = BAM_TO_FASTQ.out.fastq
                    .map { f -> tuple(f.baseName.replaceAll(/\.(fastq|fq)(\.gz)?$/, ''), f) }

                fastq_rows_ch = samplesheet_ch
                    .filter { row -> row.containsKey('fastq') && row.fastq && row.fastq.trim() != '' }
                    .map    { row -> tuple(row.sample, file(row.fastq, checkIfExists: true)) }

                input_reads_ch = fastq_rows_ch.mix(bam_converted_ch)
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
        .map { sample_id, viral_genomes, sams, primary_sams, stats, fasta ->
            tuple(sample_id, viral_genomes, sams, primary_sams,stats, fasta) }

    //QUALIMAP(MULTI_REFERENCE_MAPPING.out.sorted_bam)
    SELECT_BEST_REFERENCE(grouped_results)

    // every sample that entered the step
    all_ids = grouped_results.map { tuple(it[0], 'seen') }
    ok_ids  = SELECT_BEST_REFERENCE.out.best_ref_name.map { tuple(it[0], 'ok') }

    failed_ids = all_ids
        .join(ok_ids, remainder: true)
        .filter { id, seen, ok -> ok == null }
        .map { id, seen, ok -> id }

    failed_ids
        .collectFile(name: 'failed_best_reference_samples.txt', newLine: true, sort: true,
                    storeDir: "${params.outdir}")

    // Get unmasked fasta
    unmask_input = SELECT_BEST_REFERENCE.out.best_sam
        .join(SELECT_BEST_REFERENCE.out.mask_fa)
    
    UNMASK_SEQUENCES(unmask_input, unmask_script_ch.first())
    annotate_input = SELECT_BEST_REFERENCE.out.best_ref_fa
        .join(UNMASK_SEQUENCES.out.fasta)
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
