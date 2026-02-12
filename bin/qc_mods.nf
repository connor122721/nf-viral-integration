#!/bin/env nextflow

// Generic quality control assessment of reads
process FASTQC {
    tag "${sample_id}"
    publishDir "${params.outdir}/00_QualityControl/fastqc", mode: 'link'
    
    container params.containerQC

    input:
        tuple val(sample_id), path(reads)

    output:
        path "*.html", emit: html
        path "*.zip", emit: zip

    script:
        """
        fastqc -t ${params.threads} ${reads}
        """
}

// QC of mapping
process QUALIMAP {
    tag "${sample_id}"
    publishDir "${params.outdir}/00_QualityControl/qualimap", mode: 'link'

    container params.containerQC

    input:
        tuple val(sample_id), path(bam)

    output:
        path "${bam.baseName}_qualimap", emit: results

    script:
        """
        qualimap bamqc \
            -bam ${bam} \
            -outdir ${bam.baseName}_qualimap \
            --java-mem-size=4G \
            -nt ${params.threads}
        """
}

// Collects all fastqc / other stats and builds report html
process MULTIQC {
    publishDir "${params.outdir}/00_QualityControl/multiqc", mode: 'copy'
    
    container params.containerQC

    input:
        path '*'

    output:
        path "multiqc_report.html", emit: report
        path "multiqc_data", emit: data

    script:
        """
        multiqc . --force
        """
}