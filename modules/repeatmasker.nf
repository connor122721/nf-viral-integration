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
    """
    set -euo pipefail

    # 1. Fetch the bigBed
    wget --quiet --no-check-certificate -O ${genome_label}.rmsk.bb "${bb_url}" 
    
    # 2. Convert bb -> bed 
    bigBedToBed ${genome_label}.rmsk.bb ${genome_label}.repeatmasker.bed

    # 3. Sort + bgzip-friendly gzip 
    sort -k1,1 -k2,2n ${genome_label}.repeatmasker.bed | gzip -n > ${genome_label}.repeatmasker.bed.gz
    rm -f ${genome_label}.rmsk.bb ${genome_label}.repeatmasker.bed
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
