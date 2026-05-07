// modules/nf-core/lima/main.nf
// ============================================================================
// LIMA — PacBio HiFi/Iso-Seq demultiplexer.
//
// This module follows the nf-core/modules convention so it stays compatible
// with `nf-core modules update lima` if you later install the nf-core tooling.
// Canonical upstream:
//      https://github.com/nf-core/modules/tree/master/modules/nf-core/lima
//
// To pull the upstream version directly (recommended once you have nf-core CLI):
//      nf-core modules install lima
//
// This vendored copy uses the same process name (LIMA), the same input
// signature  ( tuple(meta, path(reads), path(barcodes)) ), and the same
// container/conda directives, so swapping in the upstream copy is a
// drop-in replacement.
// ============================================================================

process LIMA {
    tag   "${meta.id}"
    label 'process_medium'

    conda "bioconda::lima=2.9.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/lima:2.9.0--h9ee0642_0' :
        'staphb/lima:2.13.0' }"

    input:
    tuple val(meta), path(reads), path(barcodes)

    output:
    tuple val(meta), path("*.bam")           , emit: bam
    tuple val(meta), path("*.bam.pbi")       , optional: true, emit: pbi
    tuple val(meta), path("*.consensusreadset.xml"), optional: true, emit: xml
    tuple val(meta), path("*.lima.summary")  , optional: true, emit: summary
    tuple val(meta), path("*.lima.counts")   , optional: true, emit: counts
    tuple val(meta), path("*.lima.report")   , optional: true, emit: report
    tuple val(meta), path("*.json")          , optional: true, emit: json
    path  "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args   ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    """
    lima \\
        $args \\
        --num-threads $task.cpus \\
        $reads \\
        $barcodes \\
        ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lima: \$( lima --version 2>&1 | head -1 | awk '{print \$2}' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.bam.pbi
    touch ${prefix}.lima.summary
    touch ${prefix}.lima.counts
    touch ${prefix}.lima.report

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lima: 2.9.0
    END_VERSIONS
    """
}
