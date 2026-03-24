#!/bin/env nextflow
nextflow.enable.dsl = 2

// Annotate integrations and generate a per-sample HTML report
process INTEGRATION_ANNOTATE {
    tag "${sample_id}"
    publishDir "${params.outdir}/04_final_results/annotations", mode: 'copy'

    container params.container_R

    input:
        tuple val(sample_id), path(csv), path(viral_fasta), path(unmasked_fa), path(input_sam)
        path annotate_script
        path blast_script
        path sample_report_script
        path gtf

    output:
        path("*combined.csv"), emit: csv
        path("*_report.html"), emit: sample_html
        path("*nwk"), emit: nwk, optional: true
        path("*txt"), emit: txt

    script:
        def sample_id_i = sample_id.replaceAll(/.gz$/, '').replaceAll(/.fastq$/, '')
        """
        # Run the annotate script
        Rscript ${annotate_script} \\
            ${csv} \\
            ${unmasked_fa} \\
            ${viral_fasta} \\
            ${gtf} \\
            ${sample_id_i} \\
            ${input_sam}

        # Extract reference name from viral fasta header
        reference_name=\$(head -n1 ${viral_fasta} | cut -f1 -d" " | sed 's/>//g' | rev | cut -f1 -d"." | rev)

        # Make tmp directory
        mkdir -p ${projectDir}/tmp

        # Run blast proviral script — output written to .viral.tsv
        Rscript ${blast_script} \\
            --prefix ${projectDir} \\
            --in ${unmasked_fa} \\
            --virus HIV \\
            --reference \${reference_name}* \\
            --out ${sample_id_i}.viral

        # Strip trailing /ccs read-ID suffixes
        sed 's|/ccs/[0-9]*|/ccs|' ${sample_id_i}.viral.csv > \\
            ${sample_id_i}.viral_tmp.csv

        # Sort data rows on the join keys (no header present at this stage)
        sort -t',' -k7,7 ${sample_id_i}_annotated.csv > sorted1_data.csv
        sort -t',' -k1,1 ${sample_id_i}.viral_tmp.csv > sorted2_data.csv

        # Join on column 7 of annotated and column 1 of viral
        join -t',' -1 7 -2 1 sorted1_data.csv sorted2_data.csv > ${sample_id_i}.combined.csv

        # Prepend the header row
        echo "READ,INSERT,INSERT_LEN,LEFT_FLANK,RIGHT_FLANK,HUMAN_CHECK,HUMAN_ALTS,viral_sequence,viral_seq_length,viral_orientation,viral_strand,alignment_score,ref_start,ref_end,percent_identity,integration_site,viral_region,chromosome,sample,gene_name,gene_id,site_key,clone_id,reads_at_site,clone_class,is_pcr_replicate,mean_within_pid,STRAND,GENE_MATCH_STRING,MATCH_TYPE,IPDA_INTACT,IPDA_V2_INTACT,COMPLETE_5PRIME,N_GAPS_5PRIME,N_GAPS_3PRIME,N_GAPS_TOTAL,COMPLETE_3PRIME,EPISOME_FLAG" > header
        cat header ${sample_id_i}.combined.csv > ${sample_id_i}.combined.tmp.csv
        mv ${sample_id_i}.combined.tmp.csv ${sample_id_i}.combined.csv

        # Remove intermediates
        rm header sorted1_data.csv sorted2_data.csv 

        # Generate per-sample HTML report.
        Rscript ${sample_report_script} \\
            ${sample_id_i}.combined.csv \\
            ${sample_id_i} \\
            ${sample_id_i}_report.html

        exit 0
        """
}

// Generate consolidated HTML report across all samples
process CREATE_HTML_REPORT {
    publishDir "${params.outdir}/05_report", mode: 'copy'

    container params.container_R

    input:
        path combined_csvs   // collected *combined.csv files from INTEGRATION_ANNOTATE
        path report_script

    output:
        path("*_report.html"), emit: html

    script:
        def run_label = params.run_name ?: "viral_integration_run"
        """
        # Stage combined CSVs into a dedicated subdirectory.
        # create_report.R will read from here and never write back to these files.
        mkdir -p results_for_report
        for f in ${combined_csvs}; do
            cp "\${f}" results_for_report/
        done

        # Copy reference-selection summary files if available
        ref_sel_dir="${params.outdir}/01_reference_selection"
        if [ -d "\${ref_sel_dir}" ]; then
            find "\${ref_sel_dir}" -name "*_mapping_comparison.txt" \\
                -exec cp {} results_for_report/ \\;
            find "\${ref_sel_dir}" -name "*_detailed_metrics.txt" \\
                -exec cp {} results_for_report/ \\;
        fi

        # Run the project-wide report.
        Rscript ${report_script} \\
            results_for_report \\
            ${run_label} \\
            "${run_label}"
        """
}
