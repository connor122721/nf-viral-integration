#!/bin/env nextflow

nextflow.enable.dsl = 2

// Unmask sequences (extract HIV-aligned segments)
process UNMASK_SEQUENCES {
    publishDir "${params.outdir}/unmasked_sequences", mode: 'copy'

    container params.container

    input:
        path initial_sam
        path final_masked_fa
        path unmask_script

    output:
        path "*.unmasked.fa", emit: fasta

    script:
        def prefix = initial_sam.baseName.replaceAll(/\.initial.*/, "")
        """
        # Reverse mask to get HIV-aligned segments
        python ${unmask_script} ${initial_sam} ${final_masked_fa} > ${prefix}.unmasked.fa
        """
}

// Extract flanking sequences
process EXTRACT_FLANKS {
    publishDir "${params.outdir}/flanking_sequences", mode: 'copy'

    container params.container

    input:
        path masked_fasta
        path get_flanks_script

    output:
        path "*.flanks.fa", emit: fasta

    script:
        def prefix = masked_fasta.baseName.replaceAll(/\.masked/, "")
        """
        # Extract flanking regions
        python ${get_flanks_script} ${masked_fasta} > ${prefix}.flanks.fa
        """
}

// Map flanks to host genome
process MAP_FLANKS_TO_HOST {
    publishDir "${params.outdir}/flank_host_mapping", mode: 'copy'

    container params.container

    input:
        path flanks_fasta
        path host_genome

    output:
        path "*.flanks.sam", emit: sam
        path "*.flanks.bam", emit: bam
        path "*.flanks.bam.bai", emit: bai

    script:
        def prefix = flanks_fasta.baseName.replaceAll(/\.flanks/, "")
        """
        # Map flanks to host genome
        minimap2 \\
            -t ${params.threads} \\
            -Y \\
            -p 0 \\
            -N 10000 \\
            -ax map-hifi \\
            ${host_genome} \\
            ${flanks_fasta} > ${prefix}.flanks.sam

        # Convert to BAM and sort
        samtools view \\
            -@ ${params.threads} \\
            -bS ${prefix}.flanks.sam | \\
        samtools sort \\
            -@ ${params.threads} \\
            -o ${prefix}.flanks.bam

        # Index
        samtools index ${prefix}.flanks.bam
        """
}

// Confirm alignments by mapping original reads to host
process CONFIRM_HOST_ALIGNMENTS {
    publishDir "${params.outdir}/confirmed_host_mapping", mode: 'copy'

    container params.container

    input:
        path initial_sam
        path host_genome

    output:
        path "*.human.sam", emit: sam
        path "*.human.filtered.sam", emit: filtered_sam
        path "*.human.bam", emit: bam
        path "*.human.bam.bai", emit: bai

    script:
        def prefix = initial_sam.baseName.replaceAll(/\.initial.*/, "")
        """
        # Convert initial SAM to FASTQ
        samtools fastq ${initial_sam} > ${prefix}.temp.fastq

        # Map to host genome
        minimap2 \\
            -t ${params.threads} \\
            -Y \\
            -ax map-hifi \\
            ${host_genome} \\
            ${prefix}.temp.fastq > ${prefix}.human.sam

        # Filter alignments (remove unmapped and supplementary)
        samtools view \\
            -h \\
            -S \\
            -F 0x900 \\
            ${prefix}.human.sam > ${prefix}.human.filtered.sam

        # Convert to BAM
        samtools view \\
            -@ ${params.threads} \\
            -bS ${prefix}.human.filtered.sam | \\
        samtools sort \\
            -@ ${params.threads} \\
            -o ${prefix}.human.bam

        # Index
        samtools index ${prefix}.human.bam

        # Cleanup
        rm ${prefix}.temp.fastq
        """
}

// Combine HIV integration results
process COMBINE_RESULTS {
    publishDir "${params.outdir}/final_results", mode: 'copy'

    container params.container

    input:
        path flank_sam
        path host_sam
        path unmasked_fa
        path combine_script

    output:
        path "*.integration_sites.txt", emit: sites
        path "*.integration_summary.txt", emit: summary
        path "*.combine.log", emit: log

    script:
        def prefix = flank_sam.baseName.replaceAll(/\.flanks/, "")
        """
        # Combine results using the SMRTCap combine script
        python ${combine_script} ${prefix} > ${prefix}.combine.log

        # Generate summary
        echo "=== HIV Integration Analysis Summary ===" > ${prefix}.integration_summary.txt
        echo "Sample: ${prefix}" >> ${prefix}.integration_summary.txt
        echo "" >> ${prefix}.integration_summary.txt
        
        # Count integration sites from log
        if [ -f ${prefix}.combine.log ]; then
            cat ${prefix}.combine.log >> ${prefix}.integration_summary.txt
        fi
        """
}