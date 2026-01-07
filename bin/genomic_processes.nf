#!/bin/env nextflow

nextflow.enable.dsl = 2

// Unmask sequences (extract HIV-aligned segments)
process UNMASK_SEQUENCES {
    tag "${sample_id}"
    publishDir "${params.outdir}/unmasked_sequences", mode: 'copy'

    container params.container

    input:
        tuple val(sample_id), path(initial_sam), path(final_masked_fa)
        path unmask_script

    output:
        tuple val(sample_id), path("*.unmasked.fa"), emit: fasta

    script:
        def sample_id_i = sample_id.replaceAll(/.gz$/, '').replaceAll(/.fastq$/, '')
        """
        # Reverse mask to get HIV-aligned segments
        python ${unmask_script} ${initial_sam} \\
            ${final_masked_fa} > ${sample_id_i}.unmasked.fa

        """
}

// Extract flanking sequences
process EXTRACT_FLANKS {
    tag "${sample_id}"
    publishDir "${params.outdir}/flanking_sequences", mode: 'copy'

    container params.container

    input:
        tuple val(sample_id), path(masked_fasta)
        path get_flanks_script

    output:
        tuple val(sample_id), path("*.flanks.fa"), emit: fasta

    script:
        def sample_id_i = sample_id.replaceAll(/.gz$/, '').replaceAll(/.fastq$/, '')
        """
        # Extract flanking regions
        python ${get_flanks_script} ${masked_fasta} > ${sample_id_i}.flanks.fa
        """
}

// Map flanks to host genome
process MAP_FLANKS_TO_HOST {
    tag "${sample_id}"
    publishDir "${params.outdir}/flank_host_mapping", mode: 'copy'

    container params.container

    input:
        tuple val(sample_id), path(flanks_fasta)
        path host_genome

    output:
        tuple val(sample_id), path("*.flanks.sam"), emit: sam
        tuple val(sample_id), path("*.flanks.bam"), emit: bam
        tuple val(sample_id), path("*.flanks.bam.bai"), emit: bai

    script:
        def sample_id_i = sample_id.replaceAll(/.gz$/, '').replaceAll(/.fastq$/, '')
        """
        # Map flanks to host genome
        minimap2 \\
            -t ${params.threads} \\
            -Y \\
            -p 0 \\
            -N 10000 \\
            -ax map-hifi \\
            ${host_genome} \\
            ${flanks_fasta} > ${sample_id_i}.flanks.sam

        # Convert to BAM and sort
        samtools view \\
            -@ ${params.threads} \\
            -bS ${sample_id_i}.flanks.sam | \\
        samtools sort \\
            -@ ${params.threads} \\
            -o ${sample_id_i}.flanks.bam

        # Index
        samtools index ${sample_id_i}.flanks.bam
        """
}

// Confirm alignments by mapping original reads to host
process CONFIRM_HOST_ALIGNMENTS {
    tag "${sample_id}"
    publishDir "${params.outdir}/confirmed_host_mapping", mode: 'copy'

    container params.container

    input:
        tuple val(sample_id), 
              path(initial_sam)
        path host_genome

    output:
        tuple val(sample_id), path("*.human.sam"), emit: sam
        tuple val(sample_id), path("*.human.filtered.sam"), emit: filtered_sam
        tuple val(sample_id), path("*.human.bam"), emit: bam
        tuple val(sample_id), path("*.human.bam.bai"), emit: bai

    script:
        def sample_id_i = sample_id.replaceAll(/.gz$/, '').replaceAll(/.fastq$/, '')
        """
        # Convert initial SAM to FASTQ
        samtools fastq \\
            ${sample_id_i}.1.sam > ${sample_id_i}.temp.fastq

        # Map to host genome
        minimap2 \\
            -t ${params.threads} \\
            -Y \\
            -ax map-hifi \\
            ${host_genome} \\
            ${sample_id_i}.temp.fastq > ${sample_id_i}.human.sam

        # Filter alignments (remove unmapped and supplementary)
        samtools view \\
            -h \\
            -S \\
            -F 0x900 \\
            ${sample_id_i}.human.sam > ${sample_id_i}.human.filtered.sam

        # Convert to BAM
        samtools view \\
            -@ ${params.threads} \\
            -bS ${sample_id_i}.human.filtered.sam | \\
        samtools sort \\
            -@ ${params.threads} \\
            -o ${sample_id_i}.human.bam

        # Index
        samtools index ${sample_id_i}.human.bam

        # Cleanup
        rm ${sample_id_i}*.temp.fastq
        """
}

// Combine HIV integration results
process COMBINE_RESULTS {
    tag "${sample_id}"
    publishDir "${params.outdir}/03_final_results", mode: 'copy'

    container params.container

    input:
        tuple val(sample_id), path(flank_sam), path(host_sam), path(unmasked_fa)
        path combine_script

    output:
        tuple val(sample_id), path("*.integration_sites.txt"), emit: sites
        tuple val(sample_id), path("*.integration_summary.txt"), emit: summary
        tuple val(sample_id), path("*.combine.log"), emit: log

    script:
        def sample_id_i = sample_id.replaceAll(/.gz$/, '').replaceAll(/.fastq$/, '')
        """
        # Combine results using the SMRTCap combine script
        python ${combine_script} ${sample_id_i} > ${sample_id_i}.combine.log

        # Generate summary
        echo "=== HIV Integration Analysis Summary ===" > ${sample_id_i}.integration_summary.txt
        echo "Sample: ${sample_id_i}" >> ${sample_id_i}.integration_summary.txt
        echo "" >> ${sample_id_i}.integration_summary.txt
        
        # Count integration sites from log
        if [ -f ${sample_id_i}.combine.log ]; then
            cat ${sample_id_i}.combine.log >> ${sample_id_i}.integration_summary.txt
        fi
        """
}