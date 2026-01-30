#!/bin/env nextflow
nextflow.enable.dsl = 2

// Annotate integrations
process INTEGRATION_ANNOTATE {
    tag "${sample_id}"
    publishDir "${params.outdir}/04_final_results", mode: 'copy'

    container params.container_R

    input:
        tuple val(sample_id), path(csv)
        path annotate_script

    output:
        tuple val(sample_id), path("*.csv")

    script:
        def sample_id_i = sample_id.replaceAll(/.gz$/, '').replaceAll(/.fastq$/, '')
        """
        # Run the combine script with sample_id as argument
        Rscript ${combine_script} ${sample_id_i}
        exit 0
        """
}