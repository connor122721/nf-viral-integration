#!/bin/env nextflow
nextflow.enable.dsl = 2
 
// Annotation file converter - handles GTF, GTF.gz, GFF3, GFF3.gz, GFF, GFF.gz
process GFFCONVERT {
    publishDir "${params.outdir}/genome_files", mode: 'copy'
    container 'library://connmurr243/connmurr_viral/gffread.sif:latest'

    input:
        path(annotation)

    output:
        path("host.g*"), emit: gtf

    script:
        """
        # Detect file format and handle accordingly
        FNAME="${annotation}"

        if [[ "\${FNAME}" == *.gtf.gz ]]; then
            # Compressed GTF - just decompress
            gunzip -c "\${FNAME}" > host.gtf

        elif [[ "\${FNAME}" == *.gtf ]]; then
            # Uncompressed GTF - copy directly
            cp "\${FNAME}" host.gtf

        elif [[ "\${FNAME}" == *.gff3.gz ]] || [[ "\${FNAME}" == *.gff.gz ]]; then
            # Compressed GFF3/GFF - decompress then convert
            gunzip -c "\${FNAME}" > host.gff3
            #gunzip -c "\${FNAME}" > annotation_input.gff
            #gffread annotation_input.gff -T -F -o host.gtf
            #rm annotation_input.gff

        elif [[ "\${FNAME}" == *.gff3 ]] || [[ "\${FNAME}" == *.gff ]]; then
            # Uncompressed GFF3/GFF - convert directly
            cp "\${FNAME}" > host.gff3
            # gffread "\${FNAME}" -T -F -o host.gtf

        else
            echo "ERROR: Unrecognised annotation file format: \${FNAME}" >&2
            echo "Supported formats: .gtf, .gtf.gz, .gff, .gff.gz, .gff3, .gff3.gz" >&2
            exit 1
        fi
        """
}

// Unmask sequences (extract HIV-aligned segments)
process UNMASK_SEQUENCES {
    tag "${sample_id}"
    publishDir "${params.outdir}/final_results/${sample_id}", mode: 'copy'

    container params.container

    input:
        tuple val(sample_id), path(initial_sam), path(final_masked_fa)
        path unmask_script

    output:
        tuple val(sample_id), path("*.unmasked.fa"), emit: fasta

    script:
        def sample_id_i = sample_id.replaceAll(/.gz$/, '').replaceAll(/.fastq$/, '')
        """
        # Reverse mask to get HIV-aligned segments
        python ${unmask_script} ${initial_sam} \\
            ${final_masked_fa} > ${sample_id_i}.tmp.unmasked.fa

        # Drop records with no sequence characters
        awk 'BEGIN{RS=">"; FS="\\n"}
             NR>1{
               seq="";
               for(i=2;i<=NF;i++) seq=seq \$i;
               gsub(/[ \\t\\r]/,"",seq);
               if(length(seq)>0) printf ">%s", \$0
             }' ${sample_id_i}.tmp.unmasked.fa > ${sample_id_i}.unmasked.fa
        
        rm *.tmp.unmasked.fa
        echo "Finished unmasking!"
        """
}

// Extract flanking sequences
process EXTRACT_FLANKS {
    tag "${sample_id}"
    publishDir "${params.outdir}/03_flank_host_mapping/${sample_id}", mode: 'copy'

    container params.container

    input:
        tuple val(sample_id), path(masked_fasta)
        path get_flanks_script

    output:
        tuple val(sample_id), path("*.flanks.fa"), emit: fasta

    script:
        def sample_id_i = sample_id.replaceAll(/.gz$/, '').replaceAll(/.fastq$/, '')
        """
        # Extract flanking regions
        python ${get_flanks_script} ${masked_fasta} > ${sample_id_i}.flanks.fa
        """
}
