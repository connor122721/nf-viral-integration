#!/bin/env nextflow

/*
========================================================================================
    Viral Integration Detection Pipeline
========================================================================================
    Version: 2.2.0
    By: Connor S. Murray, PhD
    Based on: SMRTCap methodology (Smith Lab, University of Louisville)

    Workflow:
      1. Select best viral reference via competitive mapping
      2. Mask viral-aligned regions in reads
      3. Identify integration sites in human genome
========================================================================================
*/
nextflow.enable.dsl = 2

// ========================================================================================
// HELP MESSAGE
// ========================================================================================

def helpMessage() {
    log.info"""
    ========================================================================================
    VIRAL INTEGRATION DETECTION PIPELINE v2.2.0
    ========================================================================================

    Usage:
        nextflow run main_v2.nf --host_genome GRCh38.fa --viral_genomes "panel/*.fa"

    Required Arguments:
        --host_genome           Path to host genome FASTA
        --viral_genomes         Glob pattern for viral reference panel (e.g., "hiv_refs/*.fa")
    
    Input Options (choose one):
        --patient_dir           Directory containing patient BAM/FASTQ files
        --patient_bam           Single patient BAM file
        --patient_fastq         Single patient FASTQ file
        [OR run in simulation mode - see below]
    
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
        --min_mapq              Minimum mapping quality [default: 20]
    
    Output:
        --outdir                Output directory [default: ./output]
    
    Profiles:
        -profile slurm    Use Slurm scheduler with Apptainer
    
    Examples:
        # Analyze patient sample
        nextflow run viral_integration_smrt.nf \\
            --patient_bam sample.bam \\
            --host_genome hg38.fa \\
            --viral_genome HIV1.fa \\
            -profile slurm
        
        # Run simulation
        nextflow run viral_integration_smrt.nf \\
            --host_genome hg38.fa \\
            --viral_genome HIV1.fa \\
            --n_integrations 20 \\
            --depth 50 \\
            -profile slurm
    ========================================================================================
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

// ========================================================================================
// PRINT PARAMETER SUMMARY
// ========================================================================================

log.info ""
log.info "========================================================================================"
log.info "HIV SMRTCap INTEGRATION DETECTION PIPELINE"
log.info "========================================================================================"
log.info "Host genome       : ${params.host_genome}"
log.info "Viral genomes     : ${params.viral_genomes ?: params.viral_genome}"
log.info "Min MAPQ          : ${params.min_mapq}"
log.info "Output directory  : ${params.outdir}"
if (!params.patient_dir && !params.patient_bam && !params.patient_fastq) {
    log.info "Mode              : SIMULATION"
    log.info "Integrations      : ${params.n_integrations}"
    log.info "Depth             : ${params.depth}x"
}
log.info "========================================================================================"
log.info ""

// ========================================================================================
//  PROCESSES
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
        def seed_arg = params.seed ? "--seed ${params.seed}" : ""
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
        path "simulated_*.bam", emit: bam
        path "simulated.ccs.bam", emit: ccs
        path "simulated_*.maf.gz", emit: maf
        path "simulated_*.ref", emit: ref
    
    script:
        """
        # Simulate HiFi reads
        pbsim \\
            --strategy wgs \\
            --method qshmm \\
            --qshmm ${params.pbsim_model} \\
            --depth ${params.depth} \\
            --genome ${chimeric_fasta} \\
            --prefix simulated \\
            --pass-num=5 \\
            --length-mean ${params.length_mean} \\
            --accuracy-mean ${params.accuracy_mean}

        # Generate circular consensus sequences (CCS)
        ccs \\
            simulated*.bam \\
            simulated.ccs.bam \\
            --minPasses 3 \\
            --maxLength 50000
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
        samtools fastq \\
            -@ ${params.threads} \\
            ${bam_file} > \\
            ${bam_file.baseName}.fastq

        pigz ${bam_file.baseName}.fastq
        """
}

// Trim adapters with seqtk
process TRIM_READS {
    tag "${reads.baseName}"
    publishDir "${params.outdir}/trimmed_reads", mode: 'link'

    container params.container

    input:
        path reads

    output:
        path "*.trimmed.fastq.gz", emit: fastq

    script:
        def trim_start = params.trim_begin ?: 0
        def trim_end = params.trim_end ?: 0
        """
        seqtk trimfq -b ${trim_start} -e ${trim_end} \\
            ${reads} > ${reads.baseName}.trimmed.fastq

        pigz ${reads.baseName}.trimmed.fastq
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
        tuple val(sample_id), path(viral_genome), path("*.sam"), path("*.stats.txt"), emit: results

    script:
        def ref_name = viral_genome.baseName
        """
        # Map reads to this viral reference
        minimap2 \\
            -t ${params.threads} \\
            -Y \\
            -p 0 \\
            -N 10000 \\
            --sam-hit-only \\
            --secondary=yes \\
            -ax map-hifi \\
            ${viral_genome} \\
            ${reads} > ${sample_id}_vs_${ref_name}.sam
        
        # Remove PCR/artificial duplicates
        picard MarkDuplicates \\
            I=${sample_id}_vs_${ref_name}.sam \\
            O=${sample_id}_vs_${ref_name}.sam \\
            M=${output_prefix}_metrics.txt \\
            REMOVE_DUPLICATES=true \\
            ASSUME_SORT_ORDER=coordinate \\
            CREATE_INDEX=true \\
            VALIDATION_STRINGENCY=LENIENT \\
            MAX_RECORDS_IN_RAM=${params.picard_max_records}

        # Generate mapping statistics
        samtools flagstat ${sample_id}_vs_${ref_name}.sam > ${sample_id}_vs_${ref_name}.stats.txt
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
        tuple val(sample_id), path("*best_reference.fa"), emit: best_ref_fa
        tuple val(sample_id), path("*best_reference.sam"), emit: best_sam
        path "*mapping_comparison.txt", emit: comparison
        path "*detailed_metrics.txt", emit: metrics

    script:
        def sam_list = sam_files.collect{ it }.join(' ')
        def stats_list = stats_files.collect{ it }.join(' ')
        def genome_list = viral_genomes.collect{ it }.join(' ')
        """
        # Select best reference and get its name
        python ${projectDir}/bin/select_best_reference.py \\
            --sam-files ${sam_list} \\
            --stats-files ${stats_list} \\
            --viral-genomes ${genome_list} \\
            --output-best ${sample_id}_best_reference.txt \\
            --output-fa ${sample_id}_best_reference.fa \\
            --output-comparison ${sample_id}_mapping_comparison.txt \\
            --output-detailed ${sample_id}_detailed_metrics.txt

        # Get the best reference name (first line of best_reference.txt)
        BEST_REF=\$(head -n1 ${sample_id}_best_reference.txt)
        
        # Copy the corresponding SAM file
        cp ${sample_id}_vs_\${BEST_REF}.sam ${sample_id}_vs_\${BEST_REF}_best_reference.sam
        """
}

// Mask viral-aligned regions from reads
process MASK_VIRAL_REGIONS {
    tag "${sample_id}"
    publishDir "${params.outdir}/02_masked_reads/${sample_id}", mode: 'link'

    container params.container

    input:
        tuple val(sample_id), val(iteration), path(sam_file)

    output:
        tuple val(sample_id), val(iteration), path("*.masked.fa"), emit: fasta

    script:
        """
        # Mask viral-aligned regions
        python ${projectDir}/bin/mask.py ${sam_file} > ${sample_id}.masked.fa
        """
}

// ========================================================================================
// WORKFLOW
// ========================================================================================

// Import processes for integration site detection
include { UNMASK_SEQUENCES } from './bin/genomic_processes.nf'
include { EXTRACT_FLANKS } from './bin/genomic_processes.nf'
include { MAP_FLANKS_TO_HOST } from './bin/genomic_processes.nf'
include { CONFIRM_HOST_ALIGNMENTS } from './bin/genomic_processes.nf'
include { COMBINE_RESULTS } from './bin/genomic_processes.nf'

workflow {
    // ==================================================================================
    // SETUP: Input channels and scripts
    // ==================================================================================
    host_genome_ch = Channel.fromPath(params.host_genome, checkIfExists: true)
    viral_genomes_ch = Channel.fromPath(params.viral_genomes, checkIfExists: true)
    viral_genomes_list = Channel.fromPath(params.viral_genomes, checkIfExists: true).collect()

    // Script paths
    script_dir = "${projectDir}/bin"
    integration_script_ch = Channel.fromPath("${script_dir}/viral_integration.py", checkIfExists: true)
    unmask_script_ch = Channel.fromPath("${script_dir}/unmask.py", checkIfExists: true)
    get_flanks_script_ch = Channel.fromPath("${script_dir}/get_flanks.py", checkIfExists: true)
    combine_script_ch = Channel.fromPath("${script_dir}/combine_viral.py", checkIfExists: true)

    // ==================================================================================
    // STEP 0: Prepare input reads (simulation or patient data)
    // ==================================================================================
    if (params.patient_dir || params.patient_bam || params.patient_fastq) {
        // PATIENT DATA MODE
        if (params.patient_dir) {
            bam_files_ch = Channel.fromPath("${params.patient_dir}/*.bam", checkIfExists: false)
            fastq_files_ch = Channel.fromPath("${params.patient_dir}/*.{fastq,fq,fastq.gz,fq.gz}", checkIfExists: false)
            BAM_TO_FASTQ(bam_files_ch)
            input_reads_ch = BAM_TO_FASTQ.out.fastq.mix(fastq_files_ch)
        } else if (params.patient_bam) {
            BAM_TO_FASTQ(Channel.fromPath(params.patient_bam, checkIfExists: true))
            input_reads_ch = BAM_TO_FASTQ.out.fastq
        } else {
            input_reads_ch = Channel.fromPath(params.patient_fastq, checkIfExists: true)
        }
    } else {
        // SIMULATION MODE
        GENERATE_INTEGRATIONS(host_genome_ch, viral_genomes_ch.first(), integration_script_ch)
        SIMULATE_READS(GENERATE_INTEGRATIONS.out.fasta)
        BAM_TO_FASTQ(SIMULATE_READS.out.ccs)
        input_reads_ch = BAM_TO_FASTQ.out.fastq
    }

    // Trim reads and add sample ID
    TRIM_READS(input_reads_ch)
    trimmed_reads_ch = TRIM_READS.out.fastq
        .map { reads -> 
            def sample_id = reads.baseName.replaceAll(/\.trimmed$/, '')
            tuple(sample_id, reads)
        }

    // ==================================================================================
    // STEP 1: Select best viral reference via competitive mapping (PER SAMPLE)
    // ==================================================================================
    
    // Map each sample to all viral references
    MULTI_REFERENCE_MAPPING(trimmed_reads_ch, viral_genomes_ch)

    // Group results by sample_id
    grouped_results = MULTI_REFERENCE_MAPPING.out.results
        .groupTuple(by: 0)  // Group by sample_id
        .map { sample_id, viral_genomes, sams, stats ->
            tuple(sample_id, viral_genomes, sams, stats)
        }

    // Select best reference for each sample
    SELECT_BEST_REFERENCE(grouped_results)

    // ==================================================================================
    // STEP 2: Mask viral regions (single pass - non-iterative)
    // ==================================================================================
    // Note: We reuse the SAM file from the best reference mapping (already done in Step 1)
    
    initial_sam_ch = SELECT_BEST_REFERENCE.out.best_sam
    
    // Mask viral regions once
    mask_input = initial_sam_ch
        .map { sample_id, sam -> tuple(sample_id, 0, sam) }
    
    MASK_VIRAL_REGIONS(mask_input)

    // Extract just the masked FASTA for integration detection
    final_masked_fa = MASK_VIRAL_REGIONS.out.fasta
        .map { sample_id, iteration, masked_fa -> tuple(sample_id, masked_fa) }

    // ==================================================================================
    // STEP 3: Identify integration sites in human genome (PER SAMPLE)
    // ==================================================================================

    // Extract viral sequences from reads
    initial_sam_with_masked = initial_sam_ch
        .join(final_masked_fa)
        .map { sample_id, sam, masked_fa -> tuple(sample_id, sam, masked_fa) }

    initial_sam_with_masked.view()
    
    UNMASK_SEQUENCES(
        initial_sam_with_masked.map { sample_id, sam, masked_fa -> sam },
        initial_sam_with_masked.map { sample_id, sam, masked_fa -> masked_fa },
        unmask_script_ch)

    // Extract flanking host sequences around masked regions
    EXTRACT_FLANKS(final_masked_fa.map { sample_id, fa -> fa }, get_flanks_script_ch)

    // Map flanks to host genome to find integration sites
    MAP_FLANKS_TO_HOST(EXTRACT_FLANKS.out.fasta, host_genome_ch.first())

    // Confirm by mapping original reads to host
    CONFIRM_HOST_ALIGNMENTS(
        initial_sam_ch.map { sample_id, sam -> sam },
        host_genome_ch.first())

    // Combine results to identify precise integration sites
    COMBINE_RESULTS(
        MAP_FLANKS_TO_HOST.out.sam,
        CONFIRM_HOST_ALIGNMENTS.out.filtered_sam,
        UNMASK_SEQUENCES.out.fasta,
        combine_script_ch)
}

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
    log.info "  02_masked_reads/          - Reads with viral regions masked"
    log.info "  final_results/            - Integration sites in human genome"
    log.info ""
    log.info "Integration Sites:"
    log.info "  ${params.outdir}/final_results/*integration_sites.txt"
    log.info "  ${params.outdir}/final_results/*integration_summary.txt"
    log.info ""
    log.info "Supporting Evidence:"
    log.info "  ${params.outdir}/flank_host_mapping/*.flanks.bam"
    log.info "  ${params.outdir}/confirmed_host_mapping/*.human.bam"
    log.info "========================================================================================"
}
