#!/usr/bin/env nextflow

/*
========================================================================================
    VIRAL INTEGRATION SIMULATION PIPELINE - Lightweight Version
========================================================================================
    Simple pipeline: Generate integrations -> Simulate reads with pbsim3
    Version: 1.0.0
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
// VALIDATE PARAMETERS
// ========================================================================================

if (!params.host_genome) {
    log.error "ERROR: --host_genome is required"
    exit 1
}

if (!params.viral_genome) {
    log.error "ERROR: --viral_genome is required"
    exit 1
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

process GENERATE_INTEGRATIONS {
    publishDir "${params.outdir}/integrations", mode: 'copy'
    
    container 'docker://python:3.10'
    
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

process SIMULATE_READS {
    publishDir "${params.outdir}/simulated_reads", mode: 'copy'
    //container 'docker://ajslee/pbsim3:3.0.5'

    input:
        path chimeric_fasta
    
    output:
        path "simulated_*.bam", emit: bam
        path "simulated.ccs.bam", emit: ccs
        path "simulated_*.maf.gz", emit: maf
        path "simulated_*.ref", emit: ref
    
    script:
        """
        module load samtools/1.22.1

        # Simulate HiFi reads
        ~/pbsim \\
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
        ~/ccs \\
            simulated*.bam \\
            simulated.ccs.bam \\
            --minPasses 3 \\
            --maxLength 50000
        """
}

// ========================================================================================
// WORKFLOW
// ========================================================================================

workflow {
    // Input channels
    host_genome_ch = Channel.fromPath(params.host_genome)
    viral_genome_ch = Channel.fromPath(params.viral_genome)
    script_ch = Channel.fromPath("$projectDir/bin/viral_integration.py")
    
    // Generate integrations
    GENERATE_INTEGRATIONS(host_genome_ch,
                          viral_genome_ch,
                          script_ch)
    
    // Simulate reads with pbsim3
    SIMULATE_READS(GENERATE_INTEGRATIONS.out.fasta)

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