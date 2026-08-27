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
//
// ---------------------------------------------------------------------------
// OUTPUT FORMAT
// ---------------------------------------------------------------------------
// lima writes the format implied by the OUTPUT filename, and only certain
// in/out combinations are legal (lima.how/get-started):
//
//      In \ Out    XML   BAM   FASTQ   FASTA
//      XML         yes   yes   yes     yes
//      BAM         yes   yes   yes     yes
//      FASTQ        -     -    yes     yes
//      FASTA        -     -     -      yes
//
// BAM output from FASTQ input is IMPOSSIBLE, so this module mirrors the input
// extension rather than always writing .bam. Every read output is therefore
// declared `optional: true` — a BAM pool populates `bam`, a FASTQ pool
// populates `fastq`/`fastqgz`, never both.
// ============================================================================

process LIMA {
    tag   "${meta.id}"
    label 'process_medium'

    conda "bioconda::lima=2.9.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/lima%3A26.2.1--h9ee0642_0' :
        'staphb/lima:2.13.0' }"

    input:
    tuple val(meta), path(reads), path(barcodes)

    output:
    tuple val(meta), path("*.bam")                 , optional: true, emit: bam
    tuple val(meta), path("*.bam.pbi")             , optional: true, emit: pbi
    tuple val(meta), path("*.fastq")               , optional: true, emit: fastq
    tuple val(meta), path("*.fastq.gz")            , optional: true, emit: fastqgz
    tuple val(meta), path("*.fasta")               , optional: true, emit: fasta
    tuple val(meta), path("*.fasta.gz")            , optional: true, emit: fastagz
    tuple val(meta), path("*.consensusreadset.xml"), optional: true, emit: xml
    tuple val(meta), path("*.lima.summary")        , optional: true, emit: summary
    tuple val(meta), path("*.lima.counts")         , optional: true, emit: counts
    tuple val(meta), path("*.lima.report")         , optional: true, emit: report
    tuple val(meta), path("*.json")                , optional: true, emit: json
    path  "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // lima writes ONE merged output unless a --split* flag is given. Downstream
    // routing needs one file per barcode pair, named <prefix>.<bcF>--<bcR>.<ext>,
    // which is what --split-named produces. Note that --hifi-preset does NOT
    // imply splitting. Pass your own --split / --split-named / --split-subdirs
    // in ext.args to override this default.
    def split_arg = args.contains('--split') ? '' : '--split-named'

    if ( "${reads}" == "${prefix}.bam"      ) error "LIMA: input and output names are the same, set task.ext.prefix"
    if ( "${reads}" == "${prefix}.fastq"    ) error "LIMA: input and output names are the same, set task.ext.prefix"
    if ( "${reads}" == "${prefix}.fastq.gz" ) error "LIMA: input and output names are the same, set task.ext.prefix"
    if ( "${reads}" == "${prefix}.fasta"    ) error "LIMA: input and output names are the same, set task.ext.prefix"
    if ( "${reads}" == "${prefix}.fasta.gz" ) error "LIMA: input and output names are the same, set task.ext.prefix"
    """
    # Mirror the input format: BAM in -> BAM out, FASTQ in -> FASTQ out.
    # (.fq/.fa are normalised to .fastq/.fasta, which is what lima expects.)
    case "${reads}" in
        *.bam)              OUT_EXT="bam"                  ;;
        *.fastq.gz|*.fq.gz) OUT_EXT="fastq.gz"             ;;
        *.fastq|*.fq)       OUT_EXT="fastq"                ;;
        *.fasta.gz|*.fa.gz) OUT_EXT="fasta.gz"             ;;
        *.fasta|*.fa)       OUT_EXT="fasta"                ;;
        *.xml)              OUT_EXT="consensusreadset.xml" ;;
        *)
            echo "ERROR: LIMA cannot infer an output format for '${reads}'." >&2
            echo "       Expected .bam, .fastq[.gz], .fasta[.gz] or .xml." >&2
            exit 1
            ;;
    esac

    lima \\
        $args $split_arg \\
        --num-threads $task.cpus \\
        $reads \\
        $barcodes \\
        ${prefix}.\$OUT_EXT

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
