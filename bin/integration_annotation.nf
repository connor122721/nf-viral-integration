#!/bin/env nextflow
nextflow.enable.dsl = 2

// ============================================================================
// Annotate integrations and generate a per-sample HTML report
// ============================================================================
process INTEGRATION_ANNOTATE {
    tag "${sample_id}"
    publishDir "${params.outdir}/final_results/${sample_id}", mode: 'copy'
    container params.container_R

    input:
        tuple val(sample_id), path(viral_fasta), path(unmasked_fa), path(host_flanks), path(input_sam)
        path clone_calling_script
        path annotate_script
        path blast_script
        path sample_report_script
        path gtf
	    path repeats

    output:
        tuple val(sample_id), path("*annotated.csv"), emit: csv_ann
        path("*mapping_comparison.txt") 
        path("logs_intermediates/")
        path("blast_output/")
        path("fastas/")

    script:
        def sample_id_i = sample_id.replaceAll(/.gz$/, '').replaceAll(/.fastq$/, '')
        def report_genome = params.report_genome ?: 't2t'
        """
        reference_name=\$(head -n1 ${viral_fasta} | cut -f1 -d" " | sed 's/>//g' | rev  | cut -f1 -d"." | rev )
        ref_name=\$(head -n1 ${viral_fasta} | cut -f1 -d" " | sed 's/>//g' )

        samtools view ${input_sam} -b > tmp_working.bam
        samtools index tmp_working.bam
        
        # Run Perl script for clone-calling
        perl ${clone_calling_script} \\
            tmp_working.bam \\
            \${ref_name} \\
            ${sample_id_i}

        # 1. BLAST viral genes -----------------------------------------------
        mkdir -p ${projectDir}/tmp

        perl ${blast_script} \\
            --prefix ${projectDir} \\
            --in ${unmasked_fa} \\
            --virus HIV \\
            --reference \${reference_name}* \\
            --out ${sample_id_i}.viral.txt
        
        # Strip trailing /ccs read-ID suffixes from the viral hits
        sed 's/\t/,/g' ${sample_id_i}.viral.txt > ${sample_id_i}.viral.csv
        sed 's|/ccs/[0-9]*|/ccs|' ${sample_id_i}.viral.csv > ${sample_id_i}.viral_tmp.csv
        sed -i 's/CCS_READ_ID/READ/g' ${sample_id_i}.viral_tmp.csv

        # 2. PROPER MERGE
        Rscript ${projectDir}/bin/merge_annotated_viral.R \\
            --annotated *_MasterOfMasterFrame.tsv \\
            --viral ${sample_id_i}.viral_tmp.csv \\
            --annot_key_col 1 \\
            --viral_key_col 1 \\
            --out ${sample_id_i}.combined.csv

        # Annotate with gtf/repeatmasker
        Rscript ${annotate_script} \\
            ${sample_id_i}.combined.csv \\
            ${gtf} \\
            ${repeats}/*bed.gz \\
            ${sample_id_i}.annotated.csv

        # Move over intermediates/log files and PNGs to respective folders
        cp ${params.outdir}/01_reference_selection/${sample_id_i}/*_mapping_comparison.txt .
        mkdir -p logs_intermediates/
        cp ${params.outdir}/01_reference_selection/${sample_id_i}/*.pbmarkdup.log ./logs_intermediates/
        cp ${params.outdir}/01_reference_selection/${sample_id_i}/*\${ref_name}*.dups.readnames.txt ./logs_intermediates/
        mv CCS_ReadIDs* logs_intermediates/
        mv *png logs_intermediates/
        mv *combined.csv logs_intermediates/
        mv *viral.txt logs_intermediates/
        mkdir -p blast_output/
        mkdir -p fastas/
        cp ${params.outdir}/01_reference_selection/${sample_id_i}/*.final.*.fa fastas/

        # Conditionally copy output files
        if ls *_matches.fa 1> /dev/null 2>&1; then
            cp *_matches.fa blast_output/
        else 
            echo "Proviral BLAST screen did not start OR was inconclusive." > blast_output/warning_log.txt
        fi

        rm tmp*
        echo "Finished!"
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