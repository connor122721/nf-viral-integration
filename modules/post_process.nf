// modules/post_process.nf
// ----------------------------------------------------------------------------
// Post-processing sub-workflow for the nf-viral-integration_t2t pipeline.
// Wraps:
//   ASSIGN_CLONAL_IDS   -> bin/assign_clonal_ids.R
//   CIRCOS_SAMPLE       -> bin/make_circos.R  (mode = sample)
//   CIRCOS_PROJECT      -> bin/make_circos.R  (mode = project)
//   SAMPLE_HTML         -> bin/3.Create_Sample_HTML.R
//   PROJECT_HTML        -> bin/4.Create_Project_HTML.R

include { REPEATMASKER_DOWNLOAD; REPEATMASKER_OVERLAP } from './repeatmasker'

process ASSIGN_CLONAL_IDS {
    tag "clonal_ids"
    publishDir "${params.outdir}/04_final_results/clonal_tracking", mode: 'copy'

    input:
    path integrations_tsv
    path samplesheet

    output:
    path "integrations_with_clonal_id.tsv", emit: long
    path "clonal_persistence_wide.tsv", emit: wide
    path "clonal_summary.tsv", emit: summary

    script:
    """
    Rscript ${projectDir}/bin/assign_clonal_ids.R \\
        --integrations ${integrations_tsv} \\
        --samplesheet ${samplesheet} \\
        --window ${params.clonal_window_bp} \\
        --out_long integrations_with_clonal_id.tsv \\
        --out_wide clonal_persistence_wide.tsv \\
        --out_summary clonal_summary.tsv
    """
}

process CIRCOS_SAMPLE {
    tag { sample_id }
    publishDir "${params.outdir}/04_final_results/${sample_id}/circos", mode: 'copy'

    input:
    tuple val(sample_id), path(integrations_tsv)

    output:
    tuple val(sample_id), path("${sample_id}.circos.png"), emit: png

    script:
    def cyto_arg = params.cytoband_t2t ? "--cytoband ${params.cytoband_t2t}" : ""
    """
    Rscript ${projectDir}/bin/make_circos.R \\
        --integrations ${integrations_tsv} \\
        --mode sample \\
        --genome ${params.report_genome} \\
        ${cyto_arg} \\
        --label "${sample_id}" \\
        --out ${sample_id}.circos.png
    """
}

process CIRCOS_PROJECT {
    tag "project_circos"
    publishDir "${params.outdir}/04_final_results/project/circos", mode: 'copy'

    input:
    path integrations_tsv

    output:
    path "project.circos.png", emit: png

    script:
    def cyto_arg = params.cytoband_t2t ? "--cytoband ${params.cytoband_t2t}" : ""
    """
    Rscript ${projectDir}/bin/make_circos.R \\
        --integrations ${integrations_tsv} \\
        --mode project \\
        --genome ${params.report_genome} \\
        ${cyto_arg} \\
        --label "Project: all samples" \\
        --out project.circos.png
    """
}

process SAMPLE_HTML {
    tag { sample_id }
    publishDir "${params.outdir}/04_final_results/${sample_id}/report", mode: 'copy'

    input:
    tuple val(sample_id), path(integrations_tsv), path(circos_png)

    output:
    tuple val(sample_id), path("${sample_id}.report.html"), optional: true, emit: single
    tuple val(sample_id), path("${sample_id}_report"), optional: true, emit: multi

    script:
    """
    Rscript ${projectDir}/bin/3.Create_Sample_HTML.R \\
        --integrations ${integrations_tsv} \\
        --circos_png ${circos_png} \\
        --sample_id ${sample_id} \\
        --html_mode ${params.html_mode} \\
        --outdir .
    """
}

process PROJECT_HTML {
    tag "project_html"
    publishDir "${params.outdir}/04_final_results/project/report", mode: 'copy'

    input:
    path integrations_tsv
    path persistence_tsv
    path summary_tsv
    path circos_png
    val  sample_reports_dir

    output:
    path "project.report.html", optional: true, emit: single
    path "project_report", optional: true, emit: multi

    script:
    def reports_arg = sample_reports_dir ? "--sample_reports ${sample_reports_dir}" : ""
    """
    Rscript ${projectDir}/bin/4.Create_Project_HTML.R \\
        --integrations ${integrations_tsv} \\
        --persistence ${persistence_tsv} \\
        --summary ${summary_tsv} \\
        --circos_png ${circos_png} \\
        ${reports_arg} \\
        --html_mode ${params.html_mode} \\
        --outdir .
    """
}

// ----------------------------------------------------------------------------
// Sub-workflow that wires the above together.
//
// Inputs:
//   integrations_per_sample : channel of tuple(sample_id, path(integrations_tsv))
//   samplesheet             : single path to the run-level sample sheet
// Optional inputs (omit by passing Channel.empty()):
//   skip_repeats : boolean  (params.skip_repeatmasker)
// ----------------------------------------------------------------------------
workflow POST_PROCESS {
    take:
        integrations_per_sample
        samplesheet

    main:
        // 1. Concatenate per-sample integration tables for the project-level
        //    clonal-ID assignment step.
        all_integrations_ch = integrations_per_sample
            .map { sid, tsv -> tsv }
            .collectFile(name: 'all_samples.integrations.tsv',
                         keepHeader: true, skip: 1)

        // 2. RepeatMasker download (T2T + HG38) — cached by REPEATMASKER_DOWNLOAD.
        rmsk_inputs = Channel.of(
            tuple('t2t',  params.repeatmasker_url_t2t),
            tuple('hg38', params.repeatmasker_url_hg38))
        REPEATMASKER_DOWNLOAD(rmsk_inputs)

        // 3. Clonal IDs across all samples
        ASSIGN_CLONAL_IDS(all_integrations_ch, samplesheet)

        // 4. Re-split clonal-annotated integrations back per sample for circos
        //    + sample HTML.  Done by reading the long table downstream; keep
        //    it simple by passing the full annotated table to CIRCOS_SAMPLE
        //    via a per-sample filter (handled inside make_circos.R if needed,
        //    or upstream by splitting on sample_id with awk).
        per_sample_annotated = ASSIGN_CLONAL_IDS.out.long
            .splitCsv(header: true, sep: '\t')
            .map { row -> tuple(row.sample_id, row) }
            .groupTuple()
            .map { sid, rows ->
                def f = file("${workDir}/per_sample_${sid}.tsv")
                f.withWriter { w ->
                    w.writeLine(rows[0].keySet().join('\t'))
                    rows.each { r -> w.writeLine(r.values().join('\t')) }
                }
                tuple(sid, f)
            }

        CIRCOS_SAMPLE(per_sample_annotated)
        CIRCOS_PROJECT(ASSIGN_CLONAL_IDS.out.long)

        // 5. HTML reports
        sample_html_in = per_sample_annotated
            .join(CIRCOS_SAMPLE.out.png)
        SAMPLE_HTML(sample_html_in)

        PROJECT_HTML(
            ASSIGN_CLONAL_IDS.out.long,
            ASSIGN_CLONAL_IDS.out.wide,
            ASSIGN_CLONAL_IDS.out.summary,
            CIRCOS_PROJECT.out.png,
            "${params.outdir}")

    emit:
        clonal_long = ASSIGN_CLONAL_IDS.out.long
        clonal_wide = ASSIGN_CLONAL_IDS.out.wide
        clonal_summary = ASSIGN_CLONAL_IDS.out.summary
        sample_reports = SAMPLE_HTML.out.single.mix(SAMPLE_HTML.out.multi)
        project_report = PROJECT_HTML.out.single.mix(PROJECT_HTML.out.multi)
}
