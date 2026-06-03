#!/bin/env nextflow
nextflow.enable.dsl = 2

// ============================================================================
// Annotate integrations and generate a per-sample HTML report
// ============================================================================

process INTEGRATION_ANNOTATE {
    tag "${sample_id}"
    publishDir "${params.outdir}/04_final_results/${sample_id}/annotations", mode: 'copy'

    container params.container_R

    input:
        tuple val(sample_id), path(csv), path(viral_fasta), path(unmasked_fa), path(input_sam)
        path annotate_script
        path blast_script
        path sample_report_script
        path gtf
	path repeats

    output:
        tuple val(sample_id), path("*combined.csv"), emit: csv
        //path("*report.html"), emit: sample_html, optional: true
        //path("*circos.png"), emit: circos, optional: true
	path("*clone*"), optional: true
        path("*nwk"), emit: nwk, optional: true
        path("*viral.txt"), emit: txt

    script:
        def sample_id_i = sample_id.replaceAll(/.gz$/, '').replaceAll(/.fastq$/, '')
        def report_genome = params.report_genome ?: 't2t'
        """
        # 1. Annotate ---------------------------------------------------------
        Rscript ${annotate_script} \\
            ${csv} \\
            ${unmasked_fa} \\
            ${viral_fasta} \\
            ${gtf} \\
            ${sample_id_i} \\
            ${input_sam} \\
            ${repeats}/*bed.gz

        # 2. BLAST viral genes -----------------------------------------------
        reference_name=\$(head -n1 ${viral_fasta} | cut -f1 -d" " | sed 's/>//g' | rev | cut -f1 -d"." | rev)
        mkdir -p ${projectDir}/tmp

        Rscript ${blast_script} \\
            --prefix ${projectDir} \\
            --in ${unmasked_fa} \\
            --virus HIV \\
            --reference \${reference_name}* \\
            --out ${sample_id_i}.viral

        # Strip trailing /ccs read-ID suffixes from the viral hits
        sed 's|/ccs/[0-9]*|/ccs|' ${sample_id_i}.viral.csv > ${sample_id_i}.viral_tmp.csv
        sed -i 's/CCS_READ_ID/READ/g' ${sample_id_i}.viral_tmp.csv

        # 3. PROPER MERGE
        Rscript ${projectDir}/bin/merge_annotated_viral.R \\
            --annotated ${sample_id_i}_annotated.csv \\
            --viral ${sample_id_i}.viral_tmp.csv \\
            --annot_key_col 7 \\
            --viral_key_col 1 \\
            --out ${sample_id_i}.combined.csv
        """
}

// ============================================================================
// Generate consolidated HTML report across all samples (OLD - 05.05.2026)
// ============================================================================
process CREATE_HTML_REPORT {
    publishDir "${params.outdir}/05_report", mode: 'copy'
    container params.container_R

    input:
        path combined_csvs // collected *combined.csv files from INTEGRATION_ANNOTATE
        path report_script

    output:
        path("*_report.html"), emit: html

    script:
        def run_label = params.run_name ?: "viral_integration_run"
        """
        # Stage combined CSVs into a dedicated subdirectory.
        mkdir -p results_for_report
        for f in ${combined_csvs}; do
            cp "\${f}" results_for_report/
        done

        # Copy reference-selection summary files if available
        ref_sel_dir="${params.outdir}/01_reference_selection"
        if [ -d "\${ref_sel_dir}" ]; then
            find "\${ref_sel_dir}" -name "*_mapping_comparison.txt" -exec cp {} results_for_report/ \\;
            find "\${ref_sel_dir}" -name "*_detailed_metrics.txt"   -exec cp {} results_for_report/ \\;
        fi

        # Run the project-wide report.
        Rscript ${report_script} \\
            results_for_report \\
            ${run_label} \\
            "${run_label}"
        """
}
