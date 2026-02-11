#!/bin/env nextflow

/*
========================================================================================
    Viral Integration Detection Pipeline
========================================================================================
    Version: 2.3.0
    By: Connor S. Murray, PhD
    Based on: SMRTCap methodology (Smith Lab, University of Louisville)

    Workflow:
      1. Select best viral reference via competitive mapping
      2. Iteratively mask viral-aligned regions until no reads left
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
    VIRAL INTEGRATION DETECTION PIPELINE v2.3.0
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
        --min_mapq              Minimum mapping quality [default: 30]
        --max_iterations        Maximum iterative mapping cycles [default: 5]
    
    Output:
        --outdir                Output directory [default: ./output]
    
    Profiles:
        -profile slurm    Use Slurm scheduler with Apptainer/Singularity
    
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
log.info "Max iterations    : ${params.max_iterations}"
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
    tag "${reads.baseName.minus('.fastq')}"
    publishDir "${params.outdir}/trim_reads", mode: 'link'

    container params.container

    input:
        path reads

    output:
        path "*.trim.fastq.gz", emit: fastq

    script:
        def trim_start = params.trim_begin ?: 0
        def trim_end = params.trim_end ?: 0
        def base = reads.name.replaceAll(/\.fastq\.gz$/, '')  // Removes .fastq.gz
        """
        seqtk trimfq -b ${trim_start} -e ${trim_end} \\
            ${reads} > ${base}.trim.fastq

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
        tuple val(sample_id), path("*.dups.readnames.txt")

    script:
        def ref_name = viral_genome.baseName
        def sample_id_i = sample_id.replaceAll(/.bam$/, '')
        """
        # Map reads to this viral reference & remove reads < 200bps aligned
        minimap2 \\
            -t ${params.threads} \\
            -m 0 \\
            -Y \\
            -ax map-pb \\
            ${viral_genome} \\
            ${reads} | \\
            samtools view -h -F 4 -e 'length(seq)>=200' -b | \\
            samtools sort -@ ${params.threads} \\
                -o ${sample_id_i}_vs_${ref_name}.sorted.bam
        
        # Convert sorted BAM to FASTQ for pbmarkdup
        samtools fastq \\
            -@ ${params.threads} \\
            ${sample_id_i}_vs_${ref_name}.sorted.bam > \\
            ${sample_id_i}_vs_${ref_name}.fastq
        
        # Mark duplicates with pbmarkdup
        pbmarkdup \\
            ${sample_id_i}_vs_${ref_name}.fastq \\
            ${sample_id_i}_vs_${ref_name}.markdup.fastq \\
            --dup-file ${sample_id_i}_vs_${ref_name}.dups.fastq \\
            --log-file ${sample_id_i}_vs_${ref_name}.pbmarkdup.log 

        # Get read names of duplicates for downstream filtering
        grep "@" ${sample_id_i}_vs_${ref_name}.dups.fastq | sed 's/@//g' | cut -f1 -d" " > \\
            ${sample_id_i}_vs_${ref_name}.dups.readnames.txt

        # Convert final BAM to SAM
        samtools view -h \\
            ${sample_id_i}_vs_${ref_name}.sorted.bam \\
            --qname-file ^${sample_id_i}_vs_${ref_name}.dups.readnames.txt \\
            -o ${sample_id_i}_vs_${ref_name}.sorted.sam
        
        # Generate mapping statistics
        samtools flagstat \\
            ${sample_id_i}_vs_${ref_name}.sorted.sam > \\
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
        tuple val(sample_id), path("*best_reference.fa"), emit: best_ref_fa
        tuple val(sample_id), path("*best_reference.sam"), emit: best_sam
        path "*mapping_comparison.txt", emit: comparison
        path "*detailed_metrics.txt", emit: metrics

    script:
        def sam_list = sam_files.collect{ it }.join(' ')
        def stats_list = stats_files.collect{ it }.join(' ')
        def genome_list = viral_genomes.collect{ it }.join(' ')
        def sample_id_i = sample_id.replaceAll(/.bam$/, '')
        """
        # Select best reference and get its name
        python ${projectDir}/bin/select_best_reference.py \\
            --sam-files ${sam_list} \\
            --stats-files ${stats_list} \\
            --viral-genomes ${genome_list} \\
            --output-best ${sample_id_i}_best_reference.txt \\
            --output-fa ${sample_id_i}_best_reference.fa \\
            --output-comparison ${sample_id_i}_mapping_comparison.txt \\
            --output-detailed ${sample_id_i}_detailed_metrics.txt

        # Get the best reference name (first line of best_reference.txt)
        BEST_REF=\$(head -n1 ${sample_id_i}_best_reference.txt)
        
        # Copy the corresponding SAM file
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
        tuple val(sample_id), path("*.1.sam"), emit: iteration_1_sam
        tuple val(sample_id), path("*.final.masked.fa"), emit: final_masked_fa
        tuple val(sample_id), path("*.penultimate.masked.fa"), emit: penultimate_masked_fa
        tuple val(sample_id), path("*.[0-9].sam"), emit: iteration_sams
        path "*_iteration_log.txt", emit: log

    script:
        def sample_id_i = sample_id.replaceAll(/.bam$/, '')
        def ref_name = viral_genome.baseName
        """
        #!/bin/bash

        PREFIX="${sample_id_i}.dedup"
        MAX_ITER=${params.max_iterations} 
        
        echo "Starting iterative viral masking for ${sample_id_i}" > \${PREFIX}_iteration_log.txt
        echo "Viral reference: ${ref_name}" >> \${PREFIX}_iteration_log.txt
        echo "Max iterations: \${MAX_ITER}" >> \${PREFIX}_iteration_log.txt
        echo "" >> \${PREFIX}_iteration_log.txt
        
        # Initialize iteration counter
        NUM=1
        
        # Copy and filter initial SAM file
        cp ${initial_sam} \${PREFIX}.initial.\${NUM}.sam
        
        # Filter: remove secondary, unmapped, and supplementary (flag 0x904)
        samtools view -h -S -F 0x904 \${PREFIX}.initial.\${NUM}.sam > \${PREFIX}.\${NUM}.sam
        
        echo "Iteration \${NUM}:" >> \${PREFIX}_iteration_log.txt
        MAPPED_COUNT=\$(samtools view -c -F 4 \${PREFIX}.\${NUM}.sam)
        echo "  Mapped reads: \${MAPPED_COUNT}" >> \${PREFIX}_iteration_log.txt
        
        # Iterative loop: continue while there are mapped viral reads
        while [[ \$(samtools view -c -F 4 \${PREFIX}.\${NUM}.sam) -ne 0 ]] && [[ \${NUM} -lt \${MAX_ITER} ]]
        do
            # Mask viral regions in current iteration
            python ${projectDir}/bin/mask.py \${PREFIX}.\${NUM}.sam > \${PREFIX}.\${NUM}.fa
            
            # Remap masked sequences back to viral reference
            minimap2 -t ${params.threads} -Y -p 0 -N 10000 -ax map-pb \\
                ${viral_genome} \${PREFIX}.\${NUM}.fa > \${PREFIX}.initial.\$((NUM+1)).sam
            
            # Pick reads for next iteration
            python ${projectDir}/bin/pick_reads.py \\
                \${PREFIX}.\${NUM}.sam \\
                \${PREFIX}.initial.\$((NUM+1)).sam \\
                \${PREFIX}.\$((NUM+1)).sam
            
            # Increment counter
            NUM=\$((NUM+1))
            
            # Log progress
            echo "" >> \${PREFIX}_iteration_log.txt
            echo "Iteration \${NUM}:" >> \${PREFIX}_iteration_log.txt
            MAPPED_COUNT=\$(samtools view -c -F 4 \${PREFIX}.\${NUM}.sam)
            echo "Mapped reads: \${MAPPED_COUNT}" >> \${PREFIX}_iteration_log.txt
            
            # Check if no more viral reads found
            if [[ \$(samtools view -c -F 4 \${PREFIX}.initial.\${NUM}.sam) -eq 0 ]]; then 
                echo "No more viral reads found." >> \${PREFIX}_iteration_log.txt
                break
            fi
        done
        
        # Create final masked FASTA from last iteration
        python ${projectDir}/bin/mask.py \${PREFIX}.\${NUM}.sam > \${PREFIX}.\${NUM}.fa
        mv \${PREFIX}.\${NUM}.fa \${PREFIX}.final.masked.fa
        
        # Create penultimate masked FASTA (for flank extraction)
        if [ \${NUM} -gt 1 ]; then
            mv \${PREFIX}.\$((NUM-1)).fa \${PREFIX}.penultimate.masked.fa
        else
            # If only 1 iteration, penultimate = final
            cp \${PREFIX}.\${NUM}.fa \${PREFIX}.penultimate.masked.fa
        fi

        # Remove intermediates
        rm *initial*sam
        
        echo "" >> \${PREFIX}_iteration_log.txt
        echo "Iterative masking complete after \${NUM} iterations" >> \${PREFIX}_iteration_log.txt
        echo "Final masked FASTA: \${PREFIX}.final.masked.fa" >> \${PREFIX}_iteration_log.txt
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
include { INTEGRATION_ANNOTATE } from './bin/integration_annotation.nf'

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
    combine_script_ch = Channel.fromPath("${script_dir}/combine_hiv_V2b.py", checkIfExists: true)
    annotate_script_ch = Channel.fromPath("${script_dir}/simple_annotate_bam_v2.R", checkIfExists: true)
    blast_script_ch = Channel.fromPath("${script_dir}/findViralGenes.pl", checkIfExists: true)

    // ==================================================================================
    // STEP 0: Prepare input reads (simulation or patient-based data)
    // ==================================================================================
    if (params.patient_dir || params.patient_bam || params.patient_fastq) {

        // PATIENT DATA MODE
        if (params.patient_dir) {
            // Scan for BAM files
            bam_files_ch = Channel.fromPath("${params.patient_dir}/*.bam", checkIfExists: false)
                .filter { it.exists() }
                .filter { !it.name.toLowerCase().contains('unassigned') }
            
            // Scan for FASTQ files
            fastq_files_ch = Channel.fromPath("${params.patient_dir}/*.{fastq,fq,fastq.gz,fq.gz}", checkIfExists: false)
                .filter { it.exists() }
                .filter { !it.name.toLowerCase().contains('unassigned') }
            
            // Convert BAMs to FASTQ
            BAM_TO_FASTQ(bam_files_ch)
            
            // Merge converted + existing FASTQs
            input_reads_ch = BAM_TO_FASTQ.out.fastq.mix(fastq_files_ch)
                .map { file -> def sample_id = file.baseName.replaceAll(/\.(fastq|fq)(\.gz)?$/, '')
                    tuple(sample_id, file)}
            
            input_reads_ch.view()

        } else if (params.patient_bam) {
            bam_ch = Channel.fromPath(params.patient_bam, checkIfExists: true)            
            BAM_TO_FASTQ(bam_ch)
            input_reads_ch = BAM_TO_FASTQ.out.fastq
        } else {
            input_reads_ch = Channel.fromPath(params.patient_fastq, checkIfExists: true)
                .map { file -> def sample_id = file.baseName.replaceAll(/\.(fastq|fq)(\.gz)?$/, '')
                    tuple(sample_id, file)}
            }
                } else {
                // SIMULATION MODE
                GENERATE_INTEGRATIONS(host_genome_ch, viral_genomes_ch.first(), integration_script_ch)
                SIMULATE_READS(GENERATE_INTEGRATIONS.out.fasta)
                BAM_TO_FASTQ(SIMULATE_READS.out.bam)
                input_reads_ch = BAM_TO_FASTQ.out.fastq
                    .map { file -> def sample_id = file.baseName.replaceAll(/\.(fastq|fq)(\.gz)?$/, '')
                        tuple(sample_id, file)}
    }

    // ==================================================================================
    // STEP 1: Select best viral reference via competitive mapping per sample
    // ==================================================================================

    // Map each sample to all viral references
    MULTI_REFERENCE_MAPPING(input_reads_ch, 
                            viral_genomes_ch)

    // Group results by sample_id
    grouped_results = MULTI_REFERENCE_MAPPING.out.results
        .groupTuple(by: 0)  // Group by sample_id
        .map { sample_id, viral_genomes, sams, stats ->
            tuple(sample_id, viral_genomes, sams, stats)}

    // Select best reference for each sample
    SELECT_BEST_REFERENCE(grouped_results)

    // ==================================================================================
    // STEP 2: Iterative viral masking until no reads left
    // ==================================================================================
    
    // Combine initial SAM with best viral reference for iterative mapping
    iterative_input = SELECT_BEST_REFERENCE.out.best_sam
        .join(SELECT_BEST_REFERENCE.out.best_ref_fa)
        .map { sample_id, sam, ref_fa -> tuple(sample_id, sam, ref_fa) }
    
    // Run iterative mapping
    ITERATIVE_MAPPING(iterative_input)

    // ==================================================================================
    // STEP 3: Identify integration sites in host genome 
    // ==================================================================================

    // Extract viral sequences from reads (unmask)
    // Pass as tuple to maintain pairing
    unmask_input = ITERATIVE_MAPPING.out.iteration_1_sam
                    .join(ITERATIVE_MAPPING.out.final_masked_fa)

    UNMASK_SEQUENCES(unmask_input, 
                     unmask_script_ch.first())

    // Extract flanking host sequences (using penultimate masked FA)
    EXTRACT_FLANKS(ITERATIVE_MAPPING.out.penultimate_masked_fa,
                   get_flanks_script_ch.first())

    // Map flanks to host genome to find integration sites
    MAP_FLANKS_TO_HOST(EXTRACT_FLANKS.out.fasta, 
                       host_genome_ch.first())

    // Confirm by mapping original reads to host
    CONFIRM_HOST_ALIGNMENTS(ITERATIVE_MAPPING.out.iteration_1_sam,
                            host_genome_ch.first())

    // Combine results to identify precise integration sites
    // Join all three outputs by sample_id
    combine_input = MAP_FLANKS_TO_HOST.out.sam
        .join(CONFIRM_HOST_ALIGNMENTS.out.filtered_sam)
        .join(UNMASK_SEQUENCES.out.fasta)
        .join(ITERATIVE_MAPPING.out.iteration_sams)
        .map { sample_id, flank_sam, host_sam, iteration_sams, unmasked_fa ->
            tuple(sample_id, flank_sam, host_sam, iteration_sams, unmasked_fa)}

    // Run combination script
    COMBINE_RESULTS(combine_input, 
                    combine_script_ch.first())

    // ==================================================================================
    // STEP 4: Annotate integrate sites in host genome 
    // ==================================================================================
    INTEGRATION_ANNOTATE(COMBINE_RESULTS.out.csv,
                         SELECT_BEST_REFERENCE.out.best_ref_fa,
                         UNMASK_SEQUENCES.out.fasta,
                         annotate_script_ch.first(),
                         blast_script_ch.first())
}

// Logging workflow
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
    log.info "  04_final_results/            - Integration sites in human genome"
    log.info ""
    log.info "Integration Sites:"
    log.info "  ${params.outdir}/04_final_results/*integration_sites.txt"
    log.info "  ${params.outdir}/04_final_results/*integration_summary.txt"
    log.info ""
    log.info "Supporting Evidence:"
    log.info "  ${params.outdir}/03_flank_host_mapping/*.flanks.bam"
    log.info "  ${params.outdir}/03_flank_host_mapping/*.human.bam"
    log.info ""
    log.info "Iteration Logs:"
    log.info "  ${params.outdir}/02_iterative_masking/*_iteration_log.txt"
    log.info "========================================================================================"
}
