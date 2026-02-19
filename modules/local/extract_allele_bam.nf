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
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // Remove HLA_ prefix and make safe for filenames
    
    def allele_ref = allele.replace("HLA-", "")
    allele_safe = allele_ref.replace('*', '_').replace(':', '_')

    meta_out = meta + [allele: allele, allele_safe: allele_safe]

    """
    # Extract reads mapping to this specific allele
    # First, check if the allele exists in the reference
    echo "Extracting reads for allele: ${allele_ref}"
    REF_NAME=\$(samtools view -H ${bam} | grep '^@SQ' | grep -i "${allele_ref}" | cut -f2 | cut -d':' -f2 | head -1)
    echo "Found reference name: \$REF_NAME"
    if [ -n "\$REF_NAME" ]; then
        echo "Extracting reads for reference: \$REF_NAME"

        # Extract reads for this allele with POLYSOLVER-style filtering:
        # 1. Extract reads mapping to this allele
        # 2. Keep only reads where BOTH mates map to the SAME allele (RNEXT="=")
        # 3. Keep only properly paired reads (flag 0x2)
        # Note: Skip fixmate since reads are already from aligned BAM with correct mate info
        samtools view -h ${bam} \$REF_NAME | \\
            awk 'BEGIN {OFS="\\t"}
                 /^@/ {print; next}  # Print header lines
                 \$7 == "=" && and(\$2, 0x2) {print}  # RNEXT="=" AND properly paired flag
            ' | \\
            samtools view -b -o ${prefix}.${allele_safe}.bam -

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

    """
}
