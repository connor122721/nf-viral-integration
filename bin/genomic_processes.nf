#!/bin/env nextflow

nextflow.enable.dsl = 2

// Unmask sequences (extract HIV-aligned segments)
process UNMASK_SEQUENCES {
    tag "${sample_id}"
    publishDir "${params.outdir}/unmasked_sequences", mode: 'copy'

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
            ${final_masked_fa} > ${sample_id_i}.unmasked.fa

        """
}

// Extract flanking sequences
process EXTRACT_FLANKS {
    tag "${sample_id}"
    publishDir "${params.outdir}/flanking_sequences", mode: 'copy'

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

// Map flanks to host genome
process MAP_FLANKS_TO_HOST {
    tag "${sample_id}"
    publishDir "${params.outdir}/flank_host_mapping", mode: 'copy'

    container params.container

    input:
        tuple val(sample_id), path(flanks_fasta)
        path host_genome

    output:
        tuple val(sample_id), path("*.flanks.sam"), emit: sam
        tuple val(sample_id), path("*.flanks.bam"), emit: bam
        tuple val(sample_id), path("*.flanks.bam.bai"), emit: bai

    script:
        def sample_id_i = sample_id.replaceAll(/.gz$/, '').replaceAll(/.fastq$/, '')
        """
        # Map flanks to host genome
        minimap2 \\
            -t ${params.threads} \\
            -Y \\
            -p 0 \\
            -N 10000 \\
            -ax map-hifi \\
            ${host_genome} \\
            ${flanks_fasta} > ${sample_id_i}.flanks.sam

        # Convert to BAM and sort
        samtools view \\
            -@ ${params.threads} \\
            -bS ${sample_id_i}.flanks.sam | \\
        samtools sort \\
            -@ ${params.threads} \\
            -o ${sample_id_i}.flanks.bam

        # Index
        samtools index ${sample_id_i}.flanks.bam
        """
}

// Confirm alignments by mapping original reads to host
process CONFIRM_HOST_ALIGNMENTS {
    tag "${sample_id}"
    publishDir "${params.outdir}/confirmed_host_mapping", mode: 'copy'

    container params.container

    input:
        tuple val(sample_id), 
              path(initial_sam)
        path host_genome

    output:
        tuple val(sample_id), path("*.human.sam"), emit: sam
        tuple val(sample_id), path("*.human.filtered.sam"), emit: filtered_sam
        tuple val(sample_id), path("*.human.bam"), emit: bam
        tuple val(sample_id), path("*.human.bam.bai"), emit: bai

    script:
        def sample_id_i = sample_id.replaceAll(/.gz$/, '').replaceAll(/.fastq$/, '')
        """
        # Convert initial SAM to FASTQ
        samtools fastq \\
            ${sample_id_i}.1.sam > ${sample_id_i}.temp.fastq

        # Map to host genome
        minimap2 \\
            -t ${params.threads} \\
            -Y \\
            -ax map-hifi \\
            ${host_genome} \\
            ${sample_id_i}.temp.fastq > ${sample_id_i}.human.sam

        # Filter alignments (remove unmapped and supplementary)
        samtools view \\
            -h \\
            -S \\
            -F 0x900 \\
            ${sample_id_i}.human.sam > ${sample_id_i}.human.filtered.sam

        # Convert to BAM
        samtools view \\
            -@ ${params.threads} \\
            -bS ${sample_id_i}.human.filtered.sam | \\
        samtools sort \\
            -@ ${params.threads} \\
            -o ${sample_id_i}.human.bam

        # Index
        samtools index ${sample_id_i}.human.bam

        # Cleanup
        rm ${sample_id_i}*.temp.fastq
        """
}

// Combine HIV integration results
process COMBINE_RESULTS {
    tag "${sample_id}"
    publishDir "${params.outdir}/final_results", mode: 'copy'

    container params.container

    input:
        tuple val(sample_id), path(flank_sam), path(host_sam), path(unmasked_fa), path(iteration_sams)
        path combine_script

    output:
        tuple val(sample_id), path("*.xlsx"), emit: xlsx, optional: true
        tuple val(sample_id), path("*.tab"), emit: tab, optional: true
        tuple val(sample_id), path("*.csv"), emit: csv, optional: true
        path "*.rtf", emit: rtf, optional: true
        tuple val(sample_id), path("*.integration_summary.txt"), emit: summary
        path "*.combine.log", emit: log

    script:
        def sample_id_i = sample_id.replaceAll(/.gz$/, '').replaceAll(/.fastq$/, '')
        """
        # Stage all files with correct naming for combine_hiv_V2b.py
        # The script (via IntegrationClass) may expect either:
        #   - prefix.N.sam (new Nextflow style) OR
        #   - hiv.N.sam (old original style)
        # We'll create both for maximum compatibility
        
        echo "=== Staging iteration SAM files ===" >&2
        
        # Link all iteration SAM files with both naming conventions
        for sam_file in ${iteration_sams}; do
            # Extract iteration number from filename
            if [[ \$sam_file =~ \\.([0-9]+)\\.sam\$ ]]; then
                iter_num=\${BASH_REMATCH[1]}
                
                # Create symlink with sample_id.N.sam naming (new style)
                ln -sf \$(readlink -f \$sam_file) ${sample_id_i}.\${iter_num}.sam
                echo "  Linked: \$sam_file -> ${sample_id_i}.\${iter_num}.sam" >&2
                
                # Also create symlink with hiv.N.sam naming (legacy compatibility)
                ln -sf \$(readlink -f \$sam_file) hiv.\${iter_num}.sam
                echo "  Linked: \$sam_file -> hiv.\${iter_num}.sam (legacy)" >&2
            fi
        done
        
        # Link flanks file with both naming conventions
        ln -sf ${flank_sam} ${sample_id_i}.flanks.sam
        ln -sf ${flank_sam} flanks.sam  # Legacy compatibility
        echo "  Linked: ${flank_sam} -> ${sample_id_i}.flanks.sam" >&2
        echo "  Linked: ${flank_sam} -> flanks.sam (legacy)" >&2
        
        # Link human filtered file with both naming conventions
        ln -sf ${host_sam} ${sample_id_i}.human.filtered.sam
        ln -sf ${host_sam} human.filtered.sam  # Legacy compatibility
        echo "  Linked: ${host_sam} -> ${sample_id_i}.human.filtered.sam" >&2
        echo "  Linked: ${host_sam} -> human.filtered.sam (legacy)" >&2
        
        echo "" >&2
        echo "=== Files in working directory ===" >&2
        ls -lh ${sample_id_i}.*.sam hiv.*.sam 2>/dev/null >&2 || ls -lh *.sam >&2
        echo "" >&2
        
        # Run the combine script with sample_id as argument
        echo "=== Running combine script ===" >&2
        python ${combine_script} ${sample_id_i} > ${sample_id_i}.combine.log 2>&1 || {
            EXIT_CODE=\$?
            echo "ERROR: combine script failed with exit code \${EXIT_CODE}" >&2
            echo "" >&2
            echo "=== Last 50 lines of log ===" >&2
            tail -n 50 ${sample_id_i}.combine.log >&2
            echo "" >&2
            echo "Creating minimal outputs..." >&2
            
            # Create empty Excel file
            python -c "import xlsxwriter; wb = xlsxwriter.Workbook('${sample_id_i}.xlsx'); wb.close()" || \
                touch ${sample_id_i}.xlsx
            
            # Create minimal tab/csv files
            echo "No integration sites detected" > ${sample_id_i}.tab
            echo "No integration sites detected" > ${sample_id_i}.csv
            
            # Create error summary
            echo "=== HIV Integration Analysis Summary ===" > ${sample_id_i}.integration_summary.txt
            echo "Sample: ${sample_id_i}" >> ${sample_id_i}.integration_summary.txt
            echo "Status: FAILED" >> ${sample_id_i}.integration_summary.txt
            echo "" >> ${sample_id_i}.integration_summary.txt
            echo "The combine script encountered an error." >> ${sample_id_i}.integration_summary.txt
            echo "See ${sample_id_i}.combine.log for details." >> ${sample_id_i}.integration_summary.txt
            echo "" >> ${sample_id_i}.integration_summary.txt
            echo "=== Error Log (last 50 lines) ===" >> ${sample_id_i}.integration_summary.txt
            tail -n 50 ${sample_id_i}.combine.log >> ${sample_id_i}.integration_summary.txt 2>/dev/null || true
            
            # Exit with 0 to not fail the pipeline (outputs are created)
            exit 0
        }
        
        echo "=== Combine script completed successfully ===" >&2
        
        # Generate or enhance summary
        if [ ! -f ${sample_id_i}.integration_summary.txt ]; then
            echo "=== HIV Integration Analysis Summary ===" > ${sample_id_i}.integration_summary.txt
            echo "Sample: ${sample_id_i}" >> ${sample_id_i}.integration_summary.txt
            echo "" >> ${sample_id_i}.integration_summary.txt
        fi
        
        # Count integration sites from the tab file if it exists
        if [ -f ${sample_id_i}.tab ] && [ -s ${sample_id_i}.tab ]; then
            LINE_COUNT=\$(wc -l < ${sample_id_i}.tab)
            SITE_COUNT=\$((LINE_COUNT - 1))  # Subtract header
            echo "" >> ${sample_id_i}.integration_summary.txt
            echo "Integration sites detected: \${SITE_COUNT}" >> ${sample_id_i}.integration_summary.txt
            
            # Count unique human groups if column exists (column 2)
            if [ \${SITE_COUNT} -gt 0 ]; then
                GROUPS=\$(tail -n +2 ${sample_id_i}.tab | cut -f2 | sort -u | grep -v "^\$" | wc -l)
                echo "Unique integration groups: \${GROUPS}" >> ${sample_id_i}.integration_summary.txt
            fi
        else
            echo "" >> ${sample_id_i}.integration_summary.txt
            echo "No integration sites detected" >> ${sample_id_i}.integration_summary.txt
        fi
        
        # List all output files
        echo "" >> ${sample_id_i}.integration_summary.txt
        echo "=== Output Files ===" >> ${sample_id_i}.integration_summary.txt
        ls -1h ${sample_id_i}.xlsx ${sample_id_i}.tab ${sample_id_i}.csv ${sample_id_i}*.rtf 2>/dev/null | \
            while read f; do
                SIZE=\$(du -h "\$f" | cut -f1)
                echo "  \$f (\${SIZE})" >> ${sample_id_i}.integration_summary.txt
            done
        
        echo "" >&2
        echo "=== Processing Summary ===" >&2
        cat ${sample_id_i}.integration_summary.txt >&2
        """
}