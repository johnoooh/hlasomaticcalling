process EXTRACT_HLA_REGION {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::samtools=1.19.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0' :
        'biocontainers/samtools:1.19.2--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${prefix}.hla_region.bam"), path("${prefix}.hla_region.bam.bai"), emit: bam
    path "${prefix}.hla_stats.txt", emit: stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    // HLA region coordinates (extended MHC region on chr6)
    // Default: hg38 coordinates covering the entire MHC region
    // Can be overridden via ext.hla_region in modules.config
    def hla_region = task.ext.hla_region ?: 'chr6:28000000-34000000'

    // Alternative: Use multiple specific regions for HLA-A, B, C
    // This can be set via ext.use_specific_genes = true
    def use_specific = task.ext.use_specific_genes ?: false

    """
    # Count original reads
    TOTAL_READS=\$(samtools view -c ${bam})
    echo "Total reads in input BAM: \$TOTAL_READS"

    if [ "${use_specific}" = "true" ]; then
        # Extract reads from specific HLA-A, B, C loci (hg38 coordinates)
        # These coordinates include flanking regions to capture all HLA-relevant reads
        echo "Extracting reads from HLA-A, HLA-B, and HLA-C regions..."

        samtools view -b -h ${bam} \\
            chr6:29900000-29950000 \\
            chr6:31200000-31400000 \\
            chr6:31200000-31300000 \\
            ${args} \\
            -o ${prefix}.hla_region.bam
    else
        # Extract entire MHC region (default, matches POLYSOLVER approach)
        echo "Extracting reads from extended MHC region: ${hla_region}"

        samtools view -b -h ${bam} \\
            ${hla_region} \\
            ${args} \\
            -o ${prefix}.hla_region.bam
    fi

    # Index the extracted BAM
    samtools index ${prefix}.hla_region.bam

    # Count HLA region reads
    HLA_READS=\$(samtools view -c ${prefix}.hla_region.bam)
    echo "Reads in HLA region: \$HLA_READS"

    # Calculate percentage
    if [ "\$TOTAL_READS" -gt 0 ]; then
        PERCENT=\$(awk "BEGIN {printf \\"%.2f\\", (\$HLA_READS/\$TOTAL_READS)*100}")
    else
        PERCENT="0.00"
    fi

    # Generate statistics file
    cat > ${prefix}.hla_stats.txt <<EOF
Sample: ${meta.id}
HLA Region: ${hla_region}
Total reads in input: \$TOTAL_READS
Reads in HLA region: \$HLA_READS
Percentage of reads in HLA region: \$PERCENT%
EOF

    cat ${prefix}.hla_stats.txt

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
END_VERSIONS
    """
}
