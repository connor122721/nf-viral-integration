#!/bin/env nextflow

/*
========================================================================================
    HIV SMRTCap Integration Detection Pipeline
========================================================================================
    Version: 2.0.0
    By: Connor S. Murray, PhD
    Based on: SMRTCap methodology (Smith Lab, University of Louisville)
========================================================================================
*/
nextflow.enable.dsl = 2

// ========================================================================================
// HELP MESSAGE
// ========================================================================================

def helpMessage() {
    log.info"""
    ========================================================================================
    HIV SMRTCap INTEGRATION DETECTION PIPELINE
    ========================================================================================
    
    Usage:
        nextflow run viral_integration_smrt.nf --host_genome hg38.fa --viral_genome HIV1.fa
    
    Required Arguments:
        --host_genome           Path to host genome FASTA
        --viral_genomes         Glob pattern for multiple viral genomes (e.g., "references/HIV*.fa")
                                OR
        --viral_genome          Single viral genome FASTA (legacy support)
    
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
    
    SMRTCap Parameters:
        --max_iterations        Maximum iterative mapping cycles [default: 50]
        --trim_begin            Bases to trim from read start [default: 0]
        --trim_end              Bases to trim from read end [default: 0]
    
    Quality Filtering:
        --min_mapq              Minimum mapping quality [default: 0]
        --remove_secondary      Remove secondary alignments [default: true]
        --remove_supplementary  Remove supplementary alignments [default: true]
    
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
log.info "Max iterations    : ${params.max_iterations}"
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
    publishDir "${params.outdir}/simulated_reads", mode: 'copy'
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
    publishDir "${params.outdir}/converted_fastq", mode: 'copy'

    container params.container

    input:
        path bam_file

    output:
        path "*.fastq", emit: fastq

    script:
        """
        samtools fastq \\
            -@ ${params.threads} \\
            ${bam_file} > \\
            ${bam_file.baseName}.fastq
        """
}

// Trim adapters with seqtk
process TRIM_READS {
    tag "${reads.baseName}"
    publishDir "${params.outdir}/trimmed_reads", mode: 'copy'

    container params.container

    input:
        path reads

    output:
        path "*_trimmed.fastq", emit: fastq

    script:
        def output_name = "${reads.baseName}_trimmed.fastq"
        """
        # Uncompress if needed
        if [[ ${reads} == *.gz ]]; then
            zcat ${reads} | seqtk trimfq -b ${params.trim_begin} -e ${params.trim_end} - > ${output_name}
        else
            seqtk trimfq -b ${params.trim_begin} -e ${params.trim_end} ${reads} > ${output_name}
        fi
        """
}

// Initial alignment to all viral references for best genome selection
process MULTI_REFERENCE_MAPPING {
    tag "${reads.baseName}_vs_${viral_genome.baseName}"
    publishDir "${params.outdir}/multi_reference_mapping", mode: 'copy'

    container params.container

    input:
        path reads
        each path(viral_genome)

    output:
        tuple path(reads), path(viral_genome), path("*.sam"), path("*.stats.txt"), emit: results

    script:
        def prefix = "${reads.baseName}_vs_${viral_genome.baseName}"
        """
        # Map to this viral reference
        minimap2 \\
            -t ${params.threads} \\
            -m 0 \\
            -Y \\
            -ax map-hifi \\
            ${viral_genome} \\
            ${reads} > ${prefix}.sam

        # Filter for primary alignments only
        samtools view \\
            -h \\
            -F 0x904 \\
            ${prefix}.sam > ${prefix}.filtered.sam

        # Generate detailed statistics
        echo "=== ${viral_genome.baseName} Mapping Statistics ===" > ${prefix}.stats.txt
        samtools flagstat ${prefix}.filtered.sam >> ${prefix}.stats.txt
        echo "" >> ${prefix}.stats.txt
        echo "Mapped reads:" >> ${prefix}.stats.txt
        samtools view -c -F 4 ${prefix}.filtered.sam >> ${prefix}.stats.txt
        echo "" >> ${prefix}.stats.txt
        echo "MAPQ statistics:" >> ${prefix}.stats.txt
        samtools view ${prefix}.filtered.sam | awk '{print \$5}' | \\
            awk '{sum+=\$1; sumsq+=\$1*\$1; count++} END {
                if(count>0) {
                    mean=sum/count;
                    print "Mean_MAPQ: " mean;
                    print "Mapped_reads: " count;
                    print "Total_alignment_score: " sum;
                }
            }' >> ${prefix}.stats.txt
        
        # Move filtered sam to output
        mv ${prefix}.filtered.sam ${prefix}.sam
        """
}

// Select best viral reference genome based on mapping statistics
process SELECT_BEST_REFERENCE {
    publishDir "${params.outdir}/reference_selection", mode: 'copy'

    container params.container

    input:
        path stats_files

    output:
        path "best_reference.txt", emit: best_ref
        path "best_reference.fa", emit: best_ref_fa
        path "mapping_comparison.txt", emit: comparison

    script:
        """
        #!/usr/bin/env python3
        import re
        from pathlib import Path
        import shutil

        stats_files = """${stats_files}""".split()
        results = {}
        ref_paths = {}

        for stats_file in stats_files:
            # Parse stats file
            genome_name = Path(stats_file).stem.replace('_stats', '')
            # Extract reference name from filename pattern: reads_vs_refname
            parts = genome_name.split('_vs_')
            if len(parts) >= 2:
                ref_name = parts[-1]
            else:
                ref_name = genome_name
            
            with open(stats_file) as f:
                content = f.read()
            
            # Extract mapped reads count
            mapped_match = re.search(r'Mapped_reads:\\s*(\\d+)', content)
            mapq_match = re.search(r'Mean_MAPQ:\\s*([\\d.]+)', content)
            score_match = re.search(r'Total_alignment_score:\\s*([\\d.]+)', content)
            
            if mapped_match and mapq_match and score_match:
                mapped_reads = int(mapped_match.group(1))
                mean_mapq = float(mapq_match.group(1))
                total_score = float(score_match.group(1))
                
                # Composite score: prioritize mapped reads, then quality
                score = (mapped_reads * 1000) + mean_mapq
                
                results[ref_name] = {
                    'mapped': mapped_reads,
                    'mapq': mean_mapq,
                    'score': score,
                    'total_score': total_score
                }

        # Write comparison
        with open('mapping_comparison.txt', 'w') as f:
            f.write("Reference\\tMapped_Reads\\tMean_MAPQ\\tComposite_Score\\n")
            for ref, stats in sorted(results.items(), key=lambda x: x[1]['score'], reverse=True):
                f.write(f"{ref}\\t{stats['mapped']}\\t{stats['mapq']:.2f}\\t{stats['score']:.2f}\\n")

        # Select best
        if results:
            best = max(results.items(), key=lambda x: x[1]['score'])
            best_ref_name = best[0]
            
            with open('best_reference.txt', 'w') as f:
                f.write(best_ref_name)
            
            # Find and copy the best reference file
            # Look for the reference file in the work directory
            import glob
            possible_refs = glob.glob(f"../*/{best_ref_name}.fa") + glob.glob(f"../*/{best_ref_name}.fasta")
            if possible_refs:
                shutil.copy(possible_refs[0], 'best_reference.fa')
            else:
                # If not found, create a placeholder
                with open('best_reference.fa', 'w') as f:
                    f.write(f"# Best reference: {best_ref_name}\\n")
                    f.write("# Original file not found in work directory\\n")
            
            print(f"\\nBest reference selected: {best_ref_name}")
            print(f"  Mapped reads: {best[1]['mapped']}")
            print(f"  Mean MAPQ: {best[1]['mapq']:.2f}")
            print(f"  Composite score: {best[1]['score']:.2f}")
        else:
            with open('best_reference.txt', 'w') as f:
                f.write("NONE")
            with open('best_reference.fa', 'w') as f:
                f.write("# No suitable reference found\\n")
        """
}

// Get the actual best reference file path
process GET_BEST_REFERENCE_FILE {
    container params.container

    input:
        path best_ref_name
        path all_viral_refs

    output:
        path "selected_reference.fa", emit: reference

    script:
        """
        #!/bin/bash
        BEST_NAME=\$(cat ${best_ref_name})
        echo "Looking for reference matching: \${BEST_NAME}"
        
        # Find the matching reference file
        FOUND=0
        for ref_file in ${all_viral_refs}; do
            ref_basename=\$(basename \$ref_file .fa)
            ref_basename=\$(basename \$ref_basename .fasta)
            echo "Checking: \$ref_file (basename: \$ref_basename)"
            
            if [[ "\$ref_basename" == *"\${BEST_NAME}"* ]] || [[ "\${BEST_NAME}" == *"\$ref_basename"* ]]; then
                echo "Match found: \$ref_file"
                cp \$ref_file selected_reference.fa
                FOUND=1
                break
            fi
        done
        
        if [ \$FOUND -eq 0 ]; then
            echo "ERROR: Could not find reference file for \${BEST_NAME}"
            echo "Available files: ${all_viral_refs}"
            exit 1
        fi
        """
}

// Initial alignment to viral reference
process INITIAL_VIRAL_MAPPING {
    tag "${reads.baseName}"
    publishDir "${params.outdir}/initial_viral_mapping", mode: 'copy'

    container params.container

    input:
        path reads
        path viral_genome

    output:
        tuple val(1), path("*.initial.sam"), emit: sam
        path "*.stats.txt", emit: stats

    script:
        def prefix = reads.baseName
        """
        # Initial alignment to viral reference
        minimap2 \\
            -t ${params.threads} \\
            -m 0 \\
            -Y \\
            -ax map-hifi \\
            ${viral_genome} \\
            ${reads} > ${prefix}.initial.1.sam

        # Generate statistics
        echo "=== Initial Viral Mapping Statistics ===" > ${prefix}.initial.stats.txt
        samtools flagstat ${prefix}.initial.1.sam >> ${prefix}.initial.stats.txt
        echo "Mapped reads:" >> ${prefix}.initial.stats.txt
        samtools view -c -F 4 ${prefix}.initial.1.sam >> ${prefix}.initial.stats.txt
        """
}

// Filter alignments (remove secondary, supplementary, unmapped)
process FILTER_ALIGNMENTS {
    tag "iteration_${iteration}"
    publishDir "${params.outdir}/filtered_alignments/iteration_${iteration}", mode: 'copy'

    container params.container

    input:
        tuple val(iteration), path(sam_file)

    output:
        tuple val(iteration), path("*_filtered.sam"), emit: sam
        path "*.stats.txt", emit: stats

    script:
        def prefix = sam_file.baseName.replaceAll(/\.initial/, "")
        def flags = "260"  // unmapped + secondary
        if (params.remove_supplementary) {
            flags = "0x904"  // unmapped + secondary + supplementary
        }
        """
        # Filter alignments
        samtools view \\
            -h \\
            -S \\
            -F ${flags} \\
            ${sam_file} > ${prefix}_filtered.sam

        # Statistics
        echo "=== Iteration ${iteration} Filtered Statistics ===" > ${prefix}_filtered.stats.txt
        samtools flagstat ${prefix}_filtered.sam >> ${prefix}_filtered.stats.txt
        echo "Remaining mapped reads:" >> ${prefix}_filtered.stats.txt
        samtools view -c -F 4 ${prefix}_filtered.sam >> ${prefix}_filtered.stats.txt
        """
}

// Mask viral regions (replace with N's)
process MASK_VIRAL_REGIONS {
    tag "iteration_${iteration}"
    publishDir "${params.outdir}/masked_sequences/iteration_${iteration}", mode: 'copy'

    container params.container

    input:
        tuple val(iteration), path(sam_file)
        path mask_script

    output:
        tuple val(iteration), path("*.masked.fa"), emit: fasta
        val iteration, emit: iteration

    script:
        def prefix = sam_file.baseName.replaceAll(/_filtered/, "")
        """
        # Mask viral-aligned regions
        python ${mask_script} ${sam_file} > ${prefix}.masked.fa
        """
}

// Iterative mapping to find all viral sequences
process ITERATIVE_VIRAL_MAPPING {
    tag "iteration_${iteration}"
    publishDir "${params.outdir}/iterative_mapping/iteration_${iteration}", mode: 'copy'

    container params.container

    input:
        tuple val(iteration), path(masked_fasta)
        path viral_genome
        path pick_reads_script

    output:
        tuple val(iteration), path("*.picked.sam"), path(masked_fasta), emit: results
        path "*.stats.txt", emit: stats

    script:
        def prefix = masked_fasta.baseName.replaceAll(/\.masked/, "")
        def next_iter = iteration + 1
        """
        # Map masked sequences to viral reference
        minimap2 \\
            -t ${params.threads} \\
            -Y \\
            -p 0 \\
            -N 10000 \\
            -ax map-hifi \\
            ${viral_genome} \\
            ${masked_fasta} > ${prefix}.iter${next_iter}.initial.sam

        # Pick reads for next iteration
        # Note: This would use pick_reads.py from the original pipeline
        # For now, we'll filter the sam file
        samtools view \\
            -h \\
            -S \\
            -F 0x904 \\
            ${prefix}.iter${next_iter}.initial.sam > ${prefix}.iter${next_iter}.picked.sam

        # Statistics
        echo "=== Iteration ${next_iter} Statistics ===" > ${prefix}.iter${next_iter}.stats.txt
        samtools flagstat ${prefix}.iter${next_iter}.picked.sam >> ${prefix}.iter${next_iter}.stats.txt
        echo "New mapped reads:" >> ${prefix}.iter${next_iter}.stats.txt
        samtools view -c -F 4 ${prefix}.iter${next_iter}.picked.sam >> ${prefix}.iter${next_iter}.stats.txt
        """
}

// Check convergence
process CHECK_CONVERGENCE {
    tag "iteration_${iteration}"

    container params.container

    input:
        tuple val(iteration), path(sam_file), path(masked_fasta)

    output:
        tuple val(iteration), path(sam_file), path(masked_fasta), stdout, emit: results

    script:
        """
        # Count mapped reads
        samtools view -c -F 4 ${sam_file}
        """
}

// ========================================================================================
// WORKFLOW
// ========================================================================================

// Extract processes from external nf file
include { UNMASK_SEQUENCES as UNMASK_SEQUENCES } from './bin/genomic_processes.nf'
include { UNMASK_SEQUENCES as UNMASK_SEQUENCES_FIN } from './bin/genomic_processes.nf'
include { EXTRACT_FLANKS as EXTRACT_FLANKS } from './bin/genomic_processes.nf'
include { EXTRACT_FLANKS as EXTRACT_FLANKS_FIN } from './bin/genomic_processes.nf'
include { MAP_FLANKS_TO_HOST as MAP_FLANKS_TO_HOST } from './bin/genomic_processes.nf'  
include { MAP_FLANKS_TO_HOST as MAP_FLANKS_TO_HOST_FIN } from './bin/genomic_processes.nf' 
include { CONFIRM_HOST_ALIGNMENTS as CONFIRM_HOST_ALIGNMENTS } from './bin/genomic_processes.nf'
include { CONFIRM_HOST_ALIGNMENTS as CONFIRM_HOST_ALIGNMENTS_FIN } from './bin/genomic_processes.nf'
include { COMBINE_RESULTS as COMBINE_RESULTS } from './bin/genomic_processes.nf'
include { COMBINE_RESULTS as COMBINE_RESULTS_FIN } from './bin/genomic_processes.nf'

workflow {
    // Input channels
    host_genome_ch = Channel.fromPath(params.host_genome, checkIfExists: true)
    
    // Handle viral genome input - support both single and multiple references
    if (params.viral_genomes) {
        viral_genomes_ch = Channel.fromPath(params.viral_genomes, checkIfExists: true)
        viral_genomes_list = Channel.fromPath(params.viral_genomes, checkIfExists: true).collect()
    } else if (params.viral_genome) {
        // Legacy support for single viral genome
        viral_genomes_ch = Channel.fromPath(params.viral_genome, checkIfExists: true)
        viral_genomes_list = Channel.fromPath(params.viral_genome, checkIfExists: true).collect()
    } else {
        error "Please provide either --viral_genomes or --viral_genome"
    }
    
    // Script channels
    script_dir = "$projectDir/bin"
    integration_script_ch = Channel.fromPath("${script_dir}/viral_integration.py", checkIfExists: true)
    mask_script_ch = Channel.fromPath("${script_dir}/mask.py", checkIfExists: true)
    unmask_script_ch = Channel.fromPath("${script_dir}/unmask.py", checkIfExists: true)
    get_flanks_script_ch = Channel.fromPath("${script_dir}/get_flanks.py", checkIfExists: true)
    pick_reads_script_ch = Channel.fromPath("${script_dir}/pick_reads.py", checkIfExists: true)
    combine_script_ch = Channel.fromPath("${script_dir}/combine_viral.py", checkIfExists: true)

    // Decide whether to simulate or use patient samples
    if (params.patient_dir) {
        // Process all BAM and FASTQ files in directory
        bam_files_ch = Channel.fromPath("${params.patient_dir}/*.bam", checkIfExists: false)
        fastq_files_ch = Channel.fromPath("${params.patient_dir}/*.{fastq,fq,fastq.gz,fq.gz}", checkIfExists: false)
        
        // Convert BAMs to FASTQ
        BAM_TO_FASTQ(bam_files_ch)
        
        // Combine all FASTQs
        input_reads_ch = BAM_TO_FASTQ.out.fastq.mix(fastq_files_ch)

    } else if (params.patient_bam) {
        // Single patient BAM file
        patient_bam_ch = Channel.fromPath(params.patient_bam, checkIfExists: true)
        BAM_TO_FASTQ(patient_bam_ch)
        input_reads_ch = BAM_TO_FASTQ.out.fastq

    } else if (params.patient_fastq) {
        // Single patient FASTQ file
        input_reads_ch = Channel.fromPath(params.patient_fastq, checkIfExists: true)

    } else {
        // SIMULATION MODE - use first viral genome
        sim_viral_genome_ch = viral_genomes_ch.first()
        
        GENERATE_INTEGRATIONS(
            host_genome_ch,
            sim_viral_genome_ch,
            integration_script_ch)
        
        SIMULATE_READS(GENERATE_INTEGRATIONS.out.fasta)
        BAM_TO_FASTQ(SIMULATE_READS.out.ccs)
        input_reads_ch = BAM_TO_FASTQ.out.fastq
    }

    // Trim reads
    TRIM_READS(input_reads_ch)

    // ==================================================================================
    // STEP 1: Multi-reference mapping to select best viral genome
    // ==================================================================================
    MULTI_REFERENCE_MAPPING(
        TRIM_READS.out.fastq,
        viral_genomes_ch)

    // Collect all mapping statistics
    all_stats = MULTI_REFERENCE_MAPPING.out.results
        .map { reads, ref, sam, stats -> stats }
        .collect()

    // Select best reference based on mapping statistics
    SELECT_BEST_REFERENCE(all_stats)

    // Get the actual best reference file path
    GET_BEST_REFERENCE_FILE(
        SELECT_BEST_REFERENCE.out.best_ref,
        viral_genomes_list)

    // ==================================================================================
    // STEP 2: SMRTCap pipeline with selected best reference
    // ==================================================================================
    
    // Initial viral mapping with best reference
    INITIAL_VIRAL_MAPPING(
        TRIM_READS.out.fastq,
        GET_BEST_REFERENCE_FILE.out.reference)

    // Filter initial alignments
    FILTER_ALIGNMENTS(INITIAL_VIRAL_MAPPING.out.sam)

    // Mask viral regions
    MASK_VIRAL_REGIONS(FILTER_ALIGNMENTS.out.sam,
                       mask_script_ch)

    // Iterative mapping loop with best reference
    ITERATIVE_VIRAL_MAPPING(
        MASK_VIRAL_REGIONS.out.fasta,
        GET_BEST_REFERENCE_FILE.out.reference,
        pick_reads_script_ch)

    // Check convergence
    CHECK_CONVERGENCE(ITERATIVE_VIRAL_MAPPING.out.results)

    // Filter for converged iterations
    converged_ch = CHECK_CONVERGENCE.out.results
        .map { iter, sam, fa, count -> 
            [iter, sam, fa, count.trim().toInteger()]
        }
        .filter { iter, sam, fa, count -> 
            count == 0 || iter >= params.max_iterations 
        }

    // Get final masked fasta (from last iteration)
    final_masked_ch = converged_ch
        .map { iter, sam, fa, count -> fa }
        .first()

    // Unmask sequences
    UNMASK_SEQUENCES(
        FILTER_ALIGNMENTS.out.sam.first().map { iter, sam -> sam },
        final_masked_ch,
        unmask_script_ch)

    // Extract flanks
    EXTRACT_FLANKS(
        final_masked_ch,
        get_flanks_script_ch)

    // Map flanks to host
    MAP_FLANKS_TO_HOST(
        EXTRACT_FLANKS.out.fasta,
        host_genome_ch)

    // Confirm with original reads mapped to host
    CONFIRM_HOST_ALIGNMENTS(
        FILTER_ALIGNMENTS.out.sam.first().map { iter, sam -> sam },
        host_genome_ch)

    // Combine results
    COMBINE_RESULTS(
        MAP_FLANKS_TO_HOST.out.sam,
        CONFIRM_HOST_ALIGNMENTS.out.filtered_sam,
        UNMASK_SEQUENCES.out.fasta,
        combine_script_ch)

    // Unmask sequences
    UNMASK_SEQUENCES_FIN(
        FILTER_ALIGNMENTS.out.sam.first().map { iter, sam -> sam },
        final_masked_ch,
        unmask_script_ch)

    // Extract flanks
    EXTRACT_FLANKS_FIN(
        final_masked_ch,
        get_flanks_script_ch)

    // Map flanks to host
    MAP_FLANKS_TO_HOST_FIN(
        EXTRACT_FLANKS_FIN.out.fasta,
        host_genome_ch)

    // Confirm with original reads mapped to host
    CONFIRM_HOST_ALIGNMENTS_FIN(
        FILTER_ALIGNMENTS.out.sam.first().map { iter, sam -> sam },
        host_genome_ch)

    // Combine results
    COMBINE_RESULTS_FIN(
        MAP_FLANKS_TO_HOST_FIN.out.sam,
        CONFIRM_HOST_ALIGNMENTS_FIN.out.filtered_sam,
        UNMASK_SEQUENCES_FIN.out.fasta,
        combine_script_ch)
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
    log.info "  - Best reference:      ${params.outdir}/reference_selection/best_reference.txt"
    log.info "  - Reference comparison:${params.outdir}/reference_selection/mapping_comparison.txt"
    log.info "  - Integration sites:   ${params.outdir}/final_results/*integration_sites.txt"
    log.info "  - Summary report:      ${params.outdir}/final_results/*integration_summary.txt"
    log.info "  - Unmasked sequences:  ${params.outdir}/unmasked_sequences/*.unmasked.fa"
    log.info "  - Host alignments:     ${params.outdir}/confirmed_host_mapping/*.human.bam"
    log.info "========================================================================================"
}
