#!/usr/bin/env nextflow

/*
========================================================================================
    Viral Integration Detection Pipeline
========================================================================================
    Version: 1.0.0; By: Connor S. Murray, PhD
========================================================================================
*/
nextflow.enable.dsl = 2

// ========================================================================================
// HELP MESSAGE
// ========================================================================================

def helpMessage() {
    log.info"""
    ========================================================================================
    VIRAL INTEGRATION SIMULATION PIPELINE
    ========================================================================================
    
    Usage:
        nextflow run viral_integration.nf --host_genome hg38.fa --viral_genome HIV1.fa
    
    Required Arguments:
        --host_genome           Path to host genome FASTA
        --viral_genome          Path to viral genome FASTA
    
    Integration Options:
        --n_integrations        Number of integration sites [default: 10]
        --target_chroms         Target chromosomes (space-separated) [default: all]
        --seed                  Random seed [default: random]
    
    Read Simulation (pbsim3):
        --depth                 Sequencing depth [default: 30]
        --length_mean           Mean read length [default: 15000]
        --accuracy_mean         Mean accuracy [default: 0.99]
        --pbsim_model           pbsim3 model [default: QSHMM-RSII.model]
    
    Output:
        --outdir                Output directory [default: ./output]
    
    Examples:
        # Basic run
        nextflow run viral_integration.nf \\
            --host_genome hg38.fa \\
            --viral_genome HIV1.fa \\
            -profile slurm_apptainer
        
        # Custom parameters
        nextflow run viral_integration.nf \\
            --host_genome hg38.fa \\
            --viral_genome HIV1.fa \\
            --n_integrations 20 \\
            --depth 50 \\
            --seed 42 \\
            -profile slurm_apptainer
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
log.info "VIRAL INTEGRATION SIMULATION PIPELINE"
log.info "========================================================================================"
log.info "Host genome       : ${params.host_genome}"
log.info "Viral genome      : ${params.viral_genome}"
log.info "Integrations      : ${params.n_integrations}"
log.info "Sequencing depth  : ${params.depth}x"
log.info "Read length       : ${params.length_mean} bp"
log.info "Accuracy          : ${params.accuracy_mean}"
log.info "Output directory  : ${params.outdir}"
log.info "========================================================================================"
log.info ""

// ========================================================================================
//  PROCESSES
// ========================================================================================

// Generate in silico viral integrations
process GENERATE_INTEGRATIONS {
    publishDir "${params.outdir}/integrations", mode: 'copy'
    
    //container 'docker://python:3.10'
    container '/work/c0murr09/viral_int2.sif'
    
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
    publishDir "${params.outdir}/simulated_reads", mode: 'copy'
    //container 'docker://ajslee/pbsim3:3.0.5'
    container '/work/c0murr09/viral_int2.sif'

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
    publishDir "${params.outdir}/converted_fastq", mode: 'copy'

    container '/work/c0murr09/viral_int2.sif'

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
        
        # Fast compression
        pigz \\
            -p ${params.threads} \\
            ${bam_file.baseName}.fastq
        """
}

// Trim adapters with fastp
process FASTPLONG_TRIM {
    tag "${reads.baseName}"
    publishDir "${params.outdir}/trimmed_reads", mode: 'copy'

    container '/work/c0murr09/viral_int2.sif'

    input:
        path reads

    output:
        path "*_trimmed.fastq.gz", emit: fastq
        path "*.html", emit: html
        path "*.json", emit: json

    script:
        def output_prefix = reads.baseName.replaceAll(/\.fastq$/, "")
        def adapter_arg = params.fastp_adapters ? "--adapter_fasta ${params.fastp_adapters}" : ""
        """
        fastplong \\
            -i ${reads} \\
            -o ${output_prefix}_trimmed.fastq.gz \\
            --html ${output_prefix}_fastp.html \\
            --json ${output_prefix}_fastp.json \\
            --thread ${params.threads} \\
            --qualified_quality_phred ${params.fastp_qual} \\
            --length_required ${params.fastp_min_length} \\
            ${adapter_arg}
        """
}

// Map reads to viral genome
process MAPPING_VIRAL_GENOME {
    tag "${viral_genome.baseName}"
    publishDir "${params.outdir}/viral_mapping/${viral_genome.baseName}", mode: 'copy'

    container '/work/c0murr09/viral_int2.sif'

    input:
        path reads
        each path(viral_genome)

    output:
        path "*.sam", emit: sam
        path "*.bam", emit: bam

    script:
        def input_reads = reads.name.endsWith('.bam') ? "<(samtools fastq -@ ${params.threads} ${reads})" : reads
        def output_prefix = "${reads.baseName}_vs_${viral_genome.baseName}"
        """
        # Mapping to Viral Genome
        minimap2 \\
            -t ${params.threads} \\
            -m 0 \\
            -Y \\
            -ax map-hifi \\
            ${viral_genome} \\
            ${input_reads} > \\
            ${output_prefix}.sam

        # Convert to BAM and sort (removes secondary and unmapped alignments)
        samtools view \\
            -@ ${params.threads} \\
            -h \\
            -S \\
            -F 260,0x904 \\
            -b ${output_prefix}.sam | \\
            samtools sort \\
            -@ ${params.threads} \\
            -o ${output_prefix}.bam

        # Index bam
        samtools index ${output_prefix}.bam
        """
}

// Remove PCR duplicates with Picard
process MARK_DUPLICATES {
    tag "${bam.baseName}"
    publishDir "${params.outdir}/dedup_bams", mode: 'copy'

    container '/work/c0murr09/viral_int2.sif'

    input:
        path bam

    output:
        path "*_dedup.bam", emit: bam
        path "*_dedup.bam.bai", emit: bai
        path "*_metrics.txt", emit: metrics

    script:
        def output_prefix = bam.baseName.replaceAll(/\.bam$/, "")
        """
        picard MarkDuplicates \\
            I=${bam} \\
            O=${output_prefix}_dedup.bam \\
            M=${output_prefix}_metrics.txt \\
            REMOVE_DUPLICATES=true \\
            ASSUME_SORT_ORDER=coordinate \\
            CREATE_INDEX=true \\
            VALIDATION_STRINGENCY=LENIENT \\
            MAX_RECORDS_IN_RAM=${params.picard_max_records}

        # Index the deduplicated BAM
        samtools index ${output_prefix}_dedup.bam
        """
}

// ========================================================================================
// WORKFLOW
// ========================================================================================

workflow {
    // Input channels
    host_genome_ch = Channel.fromPath(params.host_genome)
    viral_genome_ch = Channel.fromPath(params.viral_genome)
    viral_genomes_ch = Channel.fromPath(params.viral_genomes)
    script_ch = Channel.fromPath("$projectDir/bin/viral_integration.py")

    // Decide whether to simulate or use patient samples
    if (params.patient_dir) {
        // Search for BAM and FASTQ files in the patient directory
        bam_files_ch = Channel.fromPath("${params.patient_dir}/*.bam", checkIfExists: false)
        fastq_files_ch = Channel.fromPath("${params.patient_dir}/*.fastq.gz", checkIfExists: false)

        // Convert BAMs to FASTQs
        BAM_TO_FASTQ(bam_files_ch)

        // Combine converted FASTQs with existing FASTQ files
        input_reads_ch = BAM_TO_FASTQ.out.fastq.mix(fastq_files_ch)

    } else if (params.patient_bam) {
        // Legacy: single patient BAM file - convert to FASTQ first
        patient_bam_ch = Channel.fromPath(params.patient_bam)
        BAM_TO_FASTQ(patient_bam_ch)
        input_reads_ch = BAM_TO_FASTQ.out.fastq

    } else {
        // Generate integrations
        GENERATE_INTEGRATIONS(host_genome_ch,
                              viral_genome_ch,
                              script_ch)

        // Simulate reads with pbsim3
        SIMULATE_READS(GENERATE_INTEGRATIONS.out.fasta)

        // Convert simulated BAM to FASTQ
        BAM_TO_FASTQ(SIMULATE_READS.out.ccs)
        input_reads_ch = BAM_TO_FASTQ.out.fastq
    }

    // Trim adapters with fastp
    FASTPLONG_TRIM(input_reads_ch)

    // Map to all viral genomes in parallel (each genome processed separately)
    MAPPING_VIRAL_GENOME(FASTPLONG_TRIM.out.fastq,
                         viral_genomes_ch)

    // Remove PCR duplicates from mapped BAMs
    MARK_DUPLICATES(MAPPING_VIRAL_GENOME.out.bam)

}

workflow.onComplete {
    log.info ""
    log.info "========================================================================================"
    log.info "Pipeline completed!"
    log.info "========================================================================================"
    log.info "Status:           ${workflow.success ? 'SUCCESS' : 'FAILED'}"
    log.info "Duration:         ${workflow.duration}"
    log.info "Output directory: ${params.outdir}"
    log.info ""
    log.info "Results:"
    log.info "  - Chimeric genome:     ${params.outdir}/integrations/chimeric.fa"
    log.info "  - Integration sites:   ${params.outdir}/integrations/integrations.bed"
    log.info "  - Simulated reads:     ${params.outdir}/final_reads/all_reads.fastq.gz"
    log.info "========================================================================================"
}