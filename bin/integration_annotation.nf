#!/bin/env nextflow
nextflow.enable.dsl = 2

// Annotate integrations
process INTEGRATION_ANNOTATE {
    tag "${sample_id}"
    publishDir "${params.outdir}/04_final_results", mode: 'copy'

    container params.container_R

    input:
        tuple val(sample_id), path(csv)
        tuple val(sample_id), path(fasta)
        path annotate_script
        path blast_script

    output:
        tuple val(sample_id), path("*.csv")

    script:
        def sample_id_i = sample_id.replaceAll(/.gz$/, '').replaceAll(/.fastq$/, '')
        """
        # Run the annotate script with sample_id as argument
        Rscript ${annotate_script} ${sample_id_i}
        
        # Run blast perl script
        perl ${blast_script} \\
            --in ${fasta} \\
            --reference K03455.1 \\
            --virus HIV \\
            --out ${sample_id_i}.viral.csv
        
        # Join file1 and file2 on column 1 (comma-delimited)
        join -t',' -1 9 -2 1 ${csv} ${sample_id_i}.viral.csv > \\
            ${sample_id_i}.combined.csv

        # Exit job
        exit 0
        """
}
