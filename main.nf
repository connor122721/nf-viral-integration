#!/bin/env nextflow

/*
========================================================================================
    Viral Integration Detection Pipeline 
========================================================================================
    By: Connor S. Murray, PhD
    Based on: SMRTCap methodology (Smith Lab, University of Louisville)

    Workflow:
      0. Optional lima demultiplexing of multiplexed HiFi BAMs
      1. Select best viral reference via competitive mapping
      2. Iteratively mask viral-aligned regions until no reads left
      3. Identify integration sites in human genome
      4. Annotate integration sites
      5. Assign persistent clonal IDs (chr_pos), overlap with RepeatMasker
         (T2T + HG38), render circos + slim multi-page HTML reports
========================================================================================
*/
nextflow.enable.dsl = 2

// ========================================================================================
// HELP MESSAGE
// ========================================================================================

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
        [OR run in simulation mode - see below]

    Demultiplexing (lima — modules/nf-core/lima/main.nf):
        --skip_demux            Bypass lima entirely; treat every row as already
                                demultiplexed regardless of its 'demux' value.
        --lima_args             String passed verbatim to lima.
                                Default: "--hifi-preset SYMMETRIC --min-score 80"
        --default_barcode_fasta Fallback barcode FASTA if a samplesheet row sets
                                demux=true but omits 'barcode_fasta'.

    Simulation Mode:
        --n_integrations        Number of integration sites [default: 10]
        --target_chroms         Target chromosomes (space-separated) [default: all]
        --seed                  Random seed [default: random]
        --depth                 Sequencing depth [default: 30]
        --length_mean           Mean read length [default: 15000]
        --accuracy_mean         Mean accuracy [default: 0.99]

    Masking/Mapping Parameters:
        --trim_begin            Bases to trim from read start [default: 0]
        --trim_end              Bases to trim from read end [default: 0]
        --min_mapq              Minimum mapping quality [default: 30]
        --max_iterations        Maximum iterative mapping cycles [default: 5]

    Post-processing (clonal IDs, repeats, circos, HTML reports):
        --clonal_window_bp      Window for grouping integrations into one clone [default: 10]
        --clonal_use_strand     Require same strand within a clone [default: false]
        --skip_repeatmasker     Skip RepeatMasker download + overlap [default: false]
        --report_genome         Used by circos + repeats: hg38 or t2t [default: hg38]
        --html_mode             single | multi | both [default: both]
        --legacy_project_html   Also run the old CREATE_HTML_REPORT [default: false]

    Output:
        --outdir                Output directory [default: ./output]

    Profiles:
        -profile slurm,singularity    Use Slurm scheduler with Apptainer/Singularity

    Examples:
        # Multi-sample run with mixed demuxed + multiplexed inputs
        nextflow run main.nf \\
            --samples samples.csv \\
            --host_genome hg38.fa \\
            --annotation chm13v2.0.gtf \\
            --viral_genomes "hiv_refs/*.fa" \\
            -profile slurm,singularity

        # Same, but skip lima (everything is pre-demultiplexed)
        nextflow run main.nf --samples samples.csv --skip_demux ...

        # Simulation
        nextflow run main.nf --host_genome hg38.fa --viral_genomes "hiv_refs/*.fa" \\
            --n_integrations 20 --depth 30 -profile slurm,singularity
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
log.info "Min MAPQ          : ${params.min_mapq}"
log.info "Max iterations    : ${params.max_iterations}"
log.info "Demultiplex (lima): ${params.skip_demux ? 'SKIPPED' : 'enabled (per-row demux flag)'}"
log.info "Clonal window     : ±${params.clonal_window_bp} bp"
log.info "Report genome     : ${params.report_genome}"
log.info "HTML mode         : ${params.html_mode}"
log.info "Output directory  : ${params.outdir}"
if (!params.patient_dir && !params.patient_bam && !params.patient_fastq && !params.samples) {
    log.info "Mode              : SIMULATION"
    log.info "Integrations      : ${params.n_integrations}"
    log.info "Depth             : ${params.depth}x"
}
log.info "========================================================================================"
log.info ""

// ========================================================================================
//  PROCESSES (existing — unchanged from main_2.nf)
// ========================================================================================

// Generate in silico viral integrations
process GENERATE_INTEGRATIONS {
    publishDir "${params.outdir}/integrations", mode: 'copy'
    container params.container

    input:
        path host_genome
        path viral_genome
        path script

    output:
        path "chimeric.fa", emit: fasta
        path "integrations.bed", emit: bed

    script:
        def chrom_arg = params.target_chroms ? "--target-chroms ${params.target_chroms}" : ""
        def seed_arg  = params.seed          ? "--seed ${params.seed}"                   : ""
    """
    python ${script} \\
        --host-genome ${host_genome} \\
        --viral-genome ${viral_genome} \\
        --n-integrations ${params.n_integrations} \\
        --output-fasta chimeric.fa \\
        --output-bed integrations.bed \\
        ${chrom_arg} \\
        ${seed_arg}
    """
}

// Simulate HiFi reads
process SIMULATE_READS {
    publishDir "${params.outdir}/simulated_reads", mode: 'link'
    container params.container

    input:
        path chimeric_fasta

    output:
        path "simulated_*.bam",     emit: bam
        path "simulated.ccs.bam",   emit: ccs
        path "simulated_*.maf.gz",  emit: maf
        path "simulated_*.ref",     emit: ref

    script:
    """
    pbsim \\
        --strategy wgs --method qshmm --qshmm ${params.pbsim_model} \\
        --depth ${params.depth} --genome ${chimeric_fasta} \\
        --prefix simulated --pass-num=5 \\
        --length-mean ${params.length_mean} --accuracy-mean ${params.accuracy_mean}

    ccs simulated*.bam simulated.ccs.bam --minPasses 3 --maxLength 50000
    """
}

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

// Trim adapters with seqtk
process TRIM_READS {
    tag "${reads.baseName.minus('.fastq')}"
    publishDir "${params.outdir}/trim_reads", mode: 'link'
    container params.container

    input:
        path reads

    output:
        path "*.trim.fastq.gz", emit: fastq

    script:
        def trim_start = params.trim_begin ?: 0
        def trim_end   = params.trim_end   ?: 0
        def base = reads.name.replaceAll(/\.fastq\.gz$/, '')
    """
    seqtk trimfq -b ${trim_start} -e ${trim_end} ${reads} > ${base}.trim.fastq
    pigz ${base}.trim.fastq
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

    output:
        tuple val(sample_id), path(viral_genome), path("*.sam"), path("*stats.txt"), emit: results
        tuple val(sample_id), path("*.pbmarkdup.log")
        tuple val(sample_id), path("*.readnames.txt")
        tuple val(sample_id), path("*.sorted.bam"), emit: sorted_bam

    script:
        def ref_name    = viral_genome.baseName
        def sample_id_i = sample_id.replaceAll(/.bam$/, '')
    """
    minimap2 -t ${params.threads} -m 0 -Y -ax map-pb --score-N=0 \\
        ${viral_genome} ${reads} | \\
        samtools view -h -F 4 -b | \\
        samtools sort -@ ${params.threads} -o ${sample_id_i}_vs_${ref_name}.sorted.bam

    samtools view ${sample_id_i}_vs_${ref_name}.sorted.bam | cut -f1 | sort -u > \\
        ${sample_id_i}_vs_${ref_name}.allmapped.readnames.txt

    samtools fastq -@ ${params.threads} ${sample_id_i}_vs_${ref_name}.sorted.bam > \\
        ${sample_id_i}_vs_${ref_name}.fastq

    pbmarkdup ${sample_id_i}_vs_${ref_name}.fastq \\
        ${sample_id_i}_vs_${ref_name}.markdup.fastq \\
        --dup-file ${sample_id_i}_vs_${ref_name}.dups.fastq \\
        --log-file ${sample_id_i}_vs_${ref_name}.pbmarkdup.log

    grep "@" ${sample_id_i}_vs_${ref_name}.dups.fastq | sed 's/@//g' | cut -f1 -d" " > \\
        ${sample_id_i}_vs_${ref_name}.dups.readnames.txt

    samtools view -h ${sample_id_i}_vs_${ref_name}.sorted.bam \\
        --qname-file ^${sample_id_i}_vs_${ref_name}.dups.readnames.txt \\
        -o ${sample_id_i}_vs_${ref_name}.sorted.sam

    samtools flagstat ${sample_id_i}_vs_${ref_name}.sorted.sam > \\
        ${sample_id_i}_vs_${ref_name}.stats.txt
    """
}

// Select best viral reference based on alignment metrics
process SELECT_BEST_REFERENCE {
    tag "${sample_id}"
    publishDir "${params.outdir}/01_reference_selection/${sample_id}", mode: 'link'
    container params.container

    input:
        tuple val(sample_id), path(viral_genomes), path(sam_files), path(stats_files)

    output:
        tuple val(sample_id), path("*best_reference.txt"), emit: best_ref_name
        tuple val(sample_id), path("*best_reference.fa"),  emit: best_ref_fa
        tuple val(sample_id), path("*best_reference.sam"), emit: best_sam
        path "*mapping_comparison.txt", emit: comparison
        path "*detailed_metrics.txt",   emit: metrics

    script:
        def sam_list    = sam_files.collect{ it }.join(' ')
        def stats_list  = stats_files.collect{ it }.join(' ')
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

    BEST_REF=\$(head -n1 ${sample_id_i}_best_reference.txt)
    cp ${sample_id_i}_vs_\${BEST_REF}.sam ${sample_id_i}_vs_\${BEST_REF}_best_reference.sam
    """
}

// Iterative viral mapping - loops until no viral reads remain
process ITERATIVE_MAPPING {
    tag "${sample_id}"
    publishDir "${params.outdir}/02_iterative_masking/${sample_id}", mode: 'link'
    container params.container

    input:
        tuple val(sample_id), path(initial_sam), path(viral_genome)

    output:
        tuple val(sample_id), path("*.1.sam"),                    emit: iteration_1_sam
        tuple val(sample_id), path("*.final.masked.fa"),          emit: final_masked_fa
        tuple val(sample_id), path("*.penultimate.masked.fa"),    emit: penultimate_masked_fa
        tuple val(sample_id), path("*.[0-9].sam"),                emit: iteration_sams
        path "*_iteration_log.txt",                               emit: log

    script:
        def sample_id_i = sample_id.replaceAll(/.bam$/, '')
        def ref_name    = viral_genome.baseName
    """
    #!/bin/bash
    PREFIX="${sample_id_i}.dedup"
    MAX_ITER=${params.max_iterations}

    echo "Starting iterative viral masking for ${sample_id_i}" > \${PREFIX}_iteration_log.txt
    echo "Viral reference: ${ref_name}" >> \${PREFIX}_iteration_log.txt
    echo "Max iterations: \${MAX_ITER}" >> \${PREFIX}_iteration_log.txt
    echo "" >> \${PREFIX}_iteration_log.txt

    NUM=1
    cp ${initial_sam} \${PREFIX}.initial.\${NUM}.sam
    samtools view -h -S -F 0x904 \${PREFIX}.initial.\${NUM}.sam > \${PREFIX}.\${NUM}.sam

    echo "Iteration \${NUM}:" >> \${PREFIX}_iteration_log.txt
    MAPPED_COUNT=\$(samtools view -c -F 4 \${PREFIX}.\${NUM}.sam)
    echo "  Mapped reads: \${MAPPED_COUNT}" >> \${PREFIX}_iteration_log.txt

    while [[ \$(samtools view -c -F 4 \${PREFIX}.\${NUM}.sam) -ne 0 ]] && [[ \${NUM} -lt \${MAX_ITER} ]]
    do
        python ${projectDir}/bin/mask.py \${PREFIX}.\${NUM}.sam > \${PREFIX}.\${NUM}.fa
        minimap2 -t ${params.threads} -Y -p 0 -N 10000 -ax map-pb \\
            ${viral_genome} \${PREFIX}.\${NUM}.fa > \${PREFIX}.initial.\$((NUM+1)).sam
        python ${projectDir}/bin/pick_reads.py \\
            \${PREFIX}.\${NUM}.sam \\
            \${PREFIX}.initial.\$((NUM+1)).sam \\
            \${PREFIX}.\$((NUM+1)).sam
        NUM=\$((NUM+1))
        echo "" >> \${PREFIX}_iteration_log.txt
        echo "Iteration \${NUM}:" >> \${PREFIX}_iteration_log.txt
        MAPPED_COUNT=\$(samtools view -c -F 4 \${PREFIX}.\${NUM}.sam)
        echo "Mapped reads: \${MAPPED_COUNT}" >> \${PREFIX}_iteration_log.txt
        if [[ \$(samtools view -c -F 4 \${PREFIX}.initial.\${NUM}.sam) -eq 0 ]]; then
            echo "No more viral reads found." >> \${PREFIX}_iteration_log.txt
            break
        fi
    done

    python ${projectDir}/bin/mask.py \${PREFIX}.\${NUM}.sam > \${PREFIX}.\${NUM}.fa
    mv \${PREFIX}.\${NUM}.fa \${PREFIX}.final.masked.fa

    if [ \${NUM} -gt 1 ]; then
        mv \${PREFIX}.\$((NUM-1)).fa \${PREFIX}.penultimate.masked.fa
    else
        cp \${PREFIX}.\${NUM}.fa \${PREFIX}.penultimate.masked.fa
    fi

    rm *initial*sam

    echo "" >> \${PREFIX}_iteration_log.txt
    echo "Iterative masking complete after \${NUM} iterations" >> \${PREFIX}_iteration_log.txt
    echo "Final masked FASTA: \${PREFIX}.final.masked.fa" >> \${PREFIX}_iteration_log.txt
    """
}

// ========================================================================================
//  PROCESSES (new in this refactor)
// ========================================================================================

// Rename a demuxed BAM to <sample_id>.bam so downstream BAM_TO_FASTQ + sample_id
// derivation by basename keeps working unchanged.
process RENAME_DEMUXED_BAM {
    tag "${meta.id}"
    container params.container

    input:
        tuple val(meta), path(bam)

    output:
        path "${meta.id}.bam", emit: bam

    script:
    """
    # Hard-link if same fs, otherwise copy. Either way the resulting filename
    # encodes the user-chosen sample_id from the samplesheet.
    cp -L --reflink=auto ${bam} ${meta.id}.bam || cp ${bam} ${meta.id}.bam
    """
}

// ========================================================================================
//  MODULE INCLUDES
// ========================================================================================

// QC modules (existing)
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
    host_genome_ch = Channel.fromPath(params.host_genome,   checkIfExists: true)
    viral_genomes_ch = Channel.fromPath(params.viral_genomes, checkIfExists: true)
    annotation_ch = params.annotation ? Channel.fromPath(params.annotation, checkIfExists: true) : Channel.empty()

    // Script paths
    script_dir = "${projectDir}/bin"
    integration_script_ch = Channel.fromPath("${script_dir}/viral_integration.py", checkIfExists: true)
    unmask_script_ch = Channel.fromPath("${script_dir}/unmask.py", checkIfExists: true)
    get_flanks_script_ch = Channel.fromPath("${script_dir}/get_flanks.py", checkIfExists: true)
    combine_script_ch = Channel.fromPath("${script_dir}/combine_hiv_V2b.py", checkIfExists: true)
    annotate_script_ch = Channel.fromPath("${script_dir}/1.Annotate_Flank_Bam.R", checkIfExists: true)
    blast_script_ch = Channel.fromPath("${script_dir}/2.Blast_Viral_Genes.R", checkIfExists: true)
    report_script_ch = Channel.fromPath("${script_dir}/3.Create_Sample_HTML.R", checkIfExists: true)
    project_report_script_ch = Channel.fromPath("${script_dir}/4.Create_Project_HTML.R",checkIfExists: true)

    // ==================================================================================
    // STEP 0: Prepare input reads
    //         - samplesheet: optional lima demultiplex, BAM->FASTQ, FASTQ pass-through
    //         - patient_dir / patient_bam / patient_fastq: as before
    //         - else: simulation mode
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
                    .map    { row -> tuple(row.sample, file(row.fastq, checkIfExists: true)) }

                bam_files_ch = RENAME_DEMUXED_BAM.out.bam

                BAM_TO_FASTQ(bam_files_ch)
                bam_converted_ch = BAM_TO_FASTQ.out.fastq
                    .map { f -> tuple(f.baseName.replaceAll(/\.(fastq|fq)(\.gz)?$/, ''), f) }

                input_reads_ch = fastq_rows_ch.mix(bam_converted_ch)

            } else {
                // skip_demux: previous behaviour exactly
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
    } else {
        // SIMULATION MODE -----------------------------------------------------
        GENERATE_INTEGRATIONS(host_genome_ch, viral_genomes_ch.first(), integration_script_ch)
        SIMULATE_READS(GENERATE_INTEGRATIONS.out.fasta)
        BAM_TO_FASTQ(SIMULATE_READS.out.bam)
        input_reads_ch = BAM_TO_FASTQ.out.fastq
            .map { file -> def sample_id = file.baseName.replaceAll(/\.(fastq|fq)(\.gz)?$/, '')
                tuple(sample_id, file) }
    }

    // Convert/pass-through annotation file to GTF regardless of input format
    GFFCONVERT(annotation_ch)
    gtf_ch = GFFCONVERT.out.gtf.first()

    FASTQC(input_reads_ch)

    // ==================================================================================
    // STEP 1: Select best viral reference via competitive mapping per sample
    // ==================================================================================
    MULTI_REFERENCE_MAPPING(input_reads_ch, viral_genomes_ch)

    grouped_results = MULTI_REFERENCE_MAPPING.out.results
        .groupTuple(by: 0)
        .map { sample_id, viral_genomes, sams, stats ->
            tuple(sample_id, viral_genomes, sams, stats) }

    SELECT_BEST_REFERENCE(grouped_results)

    // Run QUALIMAP only on the BEST viral-reference BAM per sample.
    //   - sorted_bam channel contains every (sample, ref) BAM produced by MULTI_REFERENCE_MAPPING
    //   - SELECT_BEST_REFERENCE.out.best_ref_name carries the chosen ref name (1 file per sample,
    //     first line of which is the ref's basename, matching ${ref}.sorted.bam naming)
    // Group all per-ref bams by sample, join with the best-ref name, then pick the matching bam.
    //best_bam_per_sample = MULTI_REFERENCE_MAPPING.out.sorted_bam
    //.groupTuple(by: 0)
    //    .join(SELECT_BEST_REFERENCE.out.best_ref_name)
    //    .map { sample_id, bams, ref_txt ->
    //        def best   = ref_txt.text.trim().tokenize('\n')[0].trim()
    //        def picked = bams.find { it.name.contains("_vs_${best}.sorted.bam") }
    //        tuple(sample_id, picked)
    //    }

    //QUALIMAP(best_bam_per_sample)

    // ==================================================================================
    // STEP 2: Iterative viral masking
    // ==================================================================================
    iterative_input = SELECT_BEST_REFERENCE.out.best_sam
        .join(SELECT_BEST_REFERENCE.out.best_ref_fa)
        .map { sample_id, sam, ref_fa -> tuple(sample_id, sam, ref_fa) }

    ITERATIVE_MAPPING(iterative_input)

    // ==================================================================================
    // STEP 3: Identify integration sites in host genome
    // ==================================================================================
    unmask_input = ITERATIVE_MAPPING.out.iteration_1_sam
        .join(ITERATIVE_MAPPING.out.final_masked_fa)
    UNMASK_SEQUENCES(unmask_input, unmask_script_ch.first())

    EXTRACT_FLANKS(ITERATIVE_MAPPING.out.penultimate_masked_fa, get_flanks_script_ch.first())

    MAP_FLANKS_TO_HOST(EXTRACT_FLANKS.out.fasta, host_genome_ch.first())

    CONFIRM_HOST_ALIGNMENTS(ITERATIVE_MAPPING.out.iteration_1_sam, host_genome_ch.first())

    combine_input = MAP_FLANKS_TO_HOST.out.sam
        .join(CONFIRM_HOST_ALIGNMENTS.out.filtered_sam)
        .join(UNMASK_SEQUENCES.out.fasta)
        .join(ITERATIVE_MAPPING.out.iteration_sams)
        .map { sample_id, flank_sam, host_sam, iteration_sams, unmasked_fa ->
            tuple(sample_id, flank_sam, host_sam, iteration_sams, unmasked_fa) }

    COMBINE_RESULTS(combine_input, combine_script_ch.first())

    // ==================================================================================
    // STEP 4: Annotate integration sites
    // ==================================================================================
    annotate_input = COMBINE_RESULTS.out.csv
        .join(SELECT_BEST_REFERENCE.out.best_ref_fa)
        .join(UNMASK_SEQUENCES.out.fasta)
        .join(ITERATIVE_MAPPING.out.iteration_1_sam)

    INTEGRATION_ANNOTATE(annotate_input,
                         annotate_script_ch.first(),
                         blast_script_ch.first(),
                         report_script_ch.first(),
                         gtf_ch)

    // Optional: keep the legacy project HTML alongside the new POST_PROCESS one.
    if (params.legacy_project_html) {
        CREATE_HTML_REPORT(INTEGRATION_ANNOTATE.out.csv.collect(),
                           project_report_script_ch.first())
    }

    // ==================================================================================
    // STEP 5: Post-processing — clonal IDs, RepeatMasker overlap, circos, slim HTML
    // ==================================================================================

    // POST_PROCESS expects per-sample integrations as tuple(sample_id, csv).
    // INTEGRATION_ANNOTATE.out.csv emits a SINGLE Path per process invocation
    // (not a tuple), so destructuring with `{ sample_id, csv -> ... }` raises
    //   Invalid method invocation `call` ... on _closure62 type
    // Derive the sample_id from the filename: "<sample>.combined.csv".
    integrations_per_sample_ch = INTEGRATION_ANNOTATE.out.csv
        .map { csv ->
            def sid = csv.getName().replaceFirst(/\.combined\.csv$/, '')
            tuple(sid, csv)
        }

    // Pass the original samplesheet so assign_clonal_ids.R can pick up
    // patient_id + timepoint columns for proper longitudinal tracking.
    samplesheet_for_post = params.samples ?
        Channel.fromPath(params.samples, checkIfExists: true) :
        Channel.empty()

    POST_PROCESS(integrations_per_sample_ch, samplesheet_for_post)

    // ==================================================================================
    // QC rollup
    // ==================================================================================
    MULTIQC(FASTQC.out.zip
        .mix(MULTI_REFERENCE_MAPPING.out.results.map { it[3] })
        .collect())
        //.mix(QUALIMAP.out.results)
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
    log.info "  02_iterative_masking/     - Iterative viral masking (until no reads left)"
    log.info "  04_final_results/         - Integration sites in human genome"
    log.info "  05_report/                - Integration sites html report/summary"
    log.info "  clonal_tracking/          - Persistent clonal IDs across timepoints"
    log.info "  <sample>/repeats/         - RepeatMasker overlap (T2T + HG38)"
    log.info "  <sample>/circos/          - Per-sample circos PNG"
    log.info "  <sample>/report/          - Slim single-page + multi-page HTML"
    log.info "  project/circos/           - Project-wide circos PNG"
    log.info "  project/report/           - Project rollup HTML (single + multi-page)"
    log.info ""
    log.info "Integration Sites:"
    log.info "  ${params.outdir}/04_final_results/*integration_sites.txt"
    log.info "  ${params.outdir}/04_final_results/*integration_summary.txt"
    log.info "  ${params.outdir}/clonal_tracking/integrations_with_clonal_id.tsv"
    log.info "  ${params.outdir}/clonal_tracking/clonal_persistence_wide.tsv"
    log.info ""
    log.info "Iteration Logs:"
    log.info "  ${params.outdir}/02_iterative_masking/*_iteration_log.txt"
    log.info "========================================================================================"
}
