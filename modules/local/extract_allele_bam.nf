process EXTRACT_ALLELE_BAM {
    tag "${meta.id}_${allele}"
    label 'process_medium'

    conda "bioconda::samtools=1.19.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0' :
        'biocontainers/samtools:1.19.2--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam), path(bai), val(allele), path(hla_reference)

    output:
    tuple val(meta_out), path("${prefix}.${allele_safe}.bam"), path("${prefix}.${allele_safe}.bam.bai"), emit: bam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // Remove HLA_ prefix and make safe for filenames
    allele_safe = allele.replaceAll("HLA_", "").replaceAll("[^A-Za-z0-9_]", "_")
    allele_ref = allele.replaceAll("HLA_", "")

    meta_out = meta + [allele: allele, allele_safe: allele_safe]

    """
    # Extract reads mapping to this specific allele
    # First, check if the allele exists in the reference
    if samtools view -H ${bam} | grep -q "${allele_ref}"; then
        REF_NAME="${allele_ref}"
    else
        # If exact match not found, try to find partial matches
        REF_NAME=\$(samtools view -H ${bam} | grep '^@SQ' | grep -i "${allele_ref}" | cut -f2 | cut -d':' -f2 | head -1)
    fi

    if [ -n "\$REF_NAME" ]; then
        echo "Extracting reads for reference: \$REF_NAME"

        # Extract reads for this allele with POLYSOLVER-style filtering:
        # 1. Extract reads mapping to this allele
        # 2. Keep only properly paired reads (both mates mapped)
        # 3. Keep only reads where BOTH mates map to the SAME allele (RNEXT="=")
        # 4. Fix mate information
        samtools view -h ${bam} \$REF_NAME | \\
            awk 'BEGIN {OFS="\\t"}
                 /^@/ {print; next}  # Print header lines
                 \$7 == "=" {print}  # Keep only reads where mate maps to same reference (RNEXT="=")
            ' | \\
            samtools view -b -h - | \\
            samtools fixmate -m - - | \\
            samtools view -b -f 0x2 - > ${prefix}.${allele_safe}.bam

        # Report statistics
        TOTAL_READS=\$(samtools view -c ${bam} \$REF_NAME || echo "0")
        FILTERED_READS=\$(samtools view -c ${prefix}.${allele_safe}.bam || echo "0")
        echo "Total reads mapping to \$REF_NAME: \$TOTAL_READS"
        echo "Filtered reads (same-allele pairs): \$FILTERED_READS"
    else
        echo "Warning: No reference found for allele ${allele_ref}"
        # Create empty BAM with header
        samtools view -H ${bam} | samtools view -b -o ${prefix}.${allele_safe}.bam -
    fi

    # Index the BAM
    samtools index ${prefix}.${allele_safe}.bam

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
END_VERSIONS
    """
}
