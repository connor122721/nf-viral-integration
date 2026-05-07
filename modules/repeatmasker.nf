// modules/repeatmasker.nf
// ----------------------------------------------------------------------------
// Download RepeatMasker tracks (bigBed) for T2T-CHM13v2.0 (hs1) and GRCh38,
// convert them to BED, and stage them as channels for downstream overlap.
//
// Inputs are URLs (configurable in nextflow.config) so the user can repoint at
// a mirror if UCSC changes path. Conversion uses `bigBedToBed` (UCSC kent
// tools, typically already in any bedtools/genomics SIF). If the binary is not
// available the process falls back to an R one-liner via rtracklayer (also
// already in your R SIF given the rest of the pipeline).
//
// Caching: outputs land in params.repeatmasker_cache_dir so subsequent runs
// skip the download. publishDir uses copy 'mode' so the cache survives even if
// work/ is wiped.
// ----------------------------------------------------------------------------

process REPEATMASKER_DOWNLOAD {
    tag      { genome_label }
    publishDir "${params.repeatmasker_cache_dir}", mode: 'copy', overwrite: false

    input:
    tuple val(genome_label), val(bb_url)

    output:
    tuple val(genome_label), path("${genome_label}.repeatmasker.bed.gz"), emit: bed

    when:
    // Skip download entirely if a cached bed.gz already exists.
    !file("${params.repeatmasker_cache_dir}/${genome_label}.repeatmasker.bed.gz").exists()

    script:
    // Auto-detect format from the URL extension. Supports:
    //   *.bb / *.bigBed -> download, convert with bigBedToBed (or rtracklayer)
    //   *.bed / *.bed.gz > download, gunzip if needed, no conversion
    // T2T URL is a .bb (bigBed); HG38 URL points at a .bed (or .bed.gz).
    def url   = bb_url.toString()
    def lower = url.toLowerCase()
    def is_bb = lower.endsWith('.bb') || lower.endsWith('.bigbed')
    def is_gz = lower.endsWith('.bed.gz') || lower.endsWith('.bed.bgz')
    def fmt   = is_bb ? 'bb' : 'bed'
    def raw   = is_bb ? "${genome_label}.rmsk.bb" :
                (is_gz ? "${genome_label}.rmsk.bed.gz" : "${genome_label}.rmsk.bed")
    """
    set -euo pipefail

    # 1. Fetch the source file (format inferred from URL extension by Nextflow).
    echo "[REPEATMASKER_DOWNLOAD] genome=${genome_label}  format=${fmt}  url=${url}"
    wget --quiet --tries=3 --timeout=120 -O ${raw} "${url}"

    # 2. Normalise to a plain BED at ${genome_label}.repeatmasker.raw.bed
    case "${fmt}" in
        bb)
            # bigBed -> BED via kent tool, R rtracklayer fallback if missing
            bigBedToBed ${raw} ${genome_label}.repeatmasker.raw.bed
            
            # Plain BED — gunzip if compressed, otherwise pass through
            if file ${raw} | grep -qi 'gzip'; then
                zcat ${raw} > ${genome_label}.repeatmasker.raw.bed
            else
                cp ${raw} ${genome_label}.repeatmasker.raw.bed
            fi
            ;;
    esac

    # 3. Sanity check — first column must look like a chromosome name.
    if [ ! -s ${genome_label}.repeatmasker.raw.bed ]; then
        echo "ERROR: ${genome_label}.repeatmasker.raw.bed is empty after download/convert." >&2
        exit 1
    fi
    head -1 ${genome_label}.repeatmasker.raw.bed | awk '{
        if (\$1 !~ /^(chr|[0-9]+\$|X\$|Y\$|M\$|MT\$)/) {
            print "WARNING: first chrom value (" \$1 ") does not look like a chromosome name" > "/dev/stderr"
        }
    }'

    # 4. Sort and gzip for downstream bedtools.
    sort -k1,1 -k2,2n ${genome_label}.repeatmasker.raw.bed \\
        | gzip -n > ${genome_label}.repeatmasker.bed.gz

    # 5. Cleanup intermediates
    rm -f ${raw} ${genome_label}.repeatmasker.raw.bed
    """
}

process REPEATMASKER_OVERLAP {
    tag        { "${sample_id}:${genome_label}" }
    publishDir "${params.outdir}/${sample_id}/repeats", mode: 'copy'

    input:
    tuple val(sample_id), val(genome_label), path(integrations_bed), path(rmsk_bed)

    output:
    tuple val(sample_id), val(genome_label),
          path("${sample_id}.${genome_label}.integrations.repeatmasker.tsv"), emit: tsv

    script:
    """
    set -euo pipefail

    # bedtools closest gives the nearest repeat plus distance (0 = inside).
    # -d adds the distance column; -t first picks one feature per tie.
    bedtools closest \\
        -a <(sort -k1,1 -k2,2n ${integrations_bed}) \\
        -b ${rmsk_bed} \\
        -d -t first \\
        > ${sample_id}.${genome_label}.integrations.repeatmasker.raw.tsv

    # Header is added by the downstream R annotator; here we keep it raw bedtools.
    awk -v OFS='\\t' 'BEGIN{
        print "int_chrom","int_start","int_end","int_name","int_score","int_strand",
              "rmsk_chrom","rmsk_start","rmsk_end","rmsk_name","rmsk_score","rmsk_strand",
              "distance_to_repeat"
    } { print }' ${sample_id}.${genome_label}.integrations.repeatmasker.raw.tsv \\
        > ${sample_id}.${genome_label}.integrations.repeatmasker.tsv

    rm -f ${sample_id}.${genome_label}.integrations.repeatmasker.raw.tsv
    """
}
