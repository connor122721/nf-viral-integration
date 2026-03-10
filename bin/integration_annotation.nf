#!/bin/env nextflow
nextflow.enable.dsl = 2

// Annotate integrations
process INTEGRATION_ANNOTATE {
    tag "${sample_id}"
    publishDir "${params.outdir}/04_final_results/annotations", mode: 'copy'

    container params.container_R

    input:
        tuple val(sample_id), path(csv)
        tuple val(sample_id), path(viral_fasta)
        tuple val(sample_id), path(unmasked_fa)
        path annotate_script
        path blast_script
        path gtf

    output:
        path("*csv"),   emit: csv
        path("*nwk"),   emit: nwk,  optional: true
        path("*txt"),   emit: txt
        path("*pdf"),   emit: pdf,  optional: true
        path("*png"),   emit: png,  optional: true

    script:
        def sample_id_i = sample_id.replaceAll(/.gz$/, '').replaceAll(/.fastq$/, '')
        """
        # Run the annotate script with sample_id as argument
        Rscript ${annotate_script} \\
            ${csv} \\
            ${unmasked_fa} \\
            ${viral_fasta} \\
            ${gtf} \\
            ${sample_id_i}
        
        # Extract reference name from viral fasta header
        reference_name=\$(head -n1 ${viral_fasta} | cut -f1 -d" " | sed 's/>//g' | rev | cut -f1 -d"." | rev)

        # Make tmp directory
        mkdir -p ${projectDir}/tmp

        # Run blast perl script
        perl ${blast_script} \\
            --prefix ${projectDir} \\
            --in ${unmasked_fa} \\
            --virus HIV \\
            --reference \${reference_name}* \\
            --out ${sample_id_i}.viral.txt

        # Convert tab-delimited viral file to csv and strip trailing /0 or /1 from read IDs
        sed 's/\t/,/g' ${sample_id_i}.viral.txt | \\
            sed 's|/ccs/[0-9]*|/ccs|' > ${sample_id_i}.viral.csv

        # Sort both files on the join column
        sort -t',' -k9,9 ${sample_id_i}_annotated.csv > sorted1.csv
        sort -t',' -k1,1 ${sample_id_i}.viral.csv > sorted2.csv

        # Join on column 9 of file1 and column 1 of file2
        join -t',' -1 9 -2 1 sorted1.csv sorted2.csv > ${sample_id_i}.combined.csv

        # Clean up
        rm sorted1.csv sorted2.csv

        # Exit job
        exit 0
        """
}

// Generate consolidated HTML report across all samples
process CREATE_HTML_REPORT {
    publishDir "${params.outdir}/05_report", mode: 'copy'

    container params.container_R

    input:
        path annotated_csvs   // collected list of *_annotated.csv files
        path report_script

    output:
        path("*_report.html"), emit: html

    script:
        def run_label = params.run_name ?: "viral_integration_run"
        """
        # Gather all annotated CSVs into a dedicated directory
        mkdir -p results_for_report
        for f in ${annotated_csvs}; do
            cp "\${f}" results_for_report/
        done

        # Generate the HTML report
        Rscript ${report_script} \\
            results_for_report \\
            ${run_label} \\
            "${run_label}"
        """
}
