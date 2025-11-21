process COMBINE_ALLELE_VCFS {
    tag "${meta.sample_id}_${caller}"
    label 'process_low'

    conda "bioconda::bcftools=1.19"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.19--h8b25389_0' :
        'biocontainers/bcftools:1.19--h8b25389_0' }"

    input:
    tuple val(meta), val(caller), path(vcfs), path(tbis)

    output:
    tuple val(meta_out), path("${prefix}.combined.vcf.gz"), path("${prefix}.combined.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.sample_id}_${caller}"
    meta_out = meta + [id: "${meta.sample_id}_${caller}_combined"]

    """
    # Create a list of VCF files to combine
    echo "${vcfs.join('\n')}" > vcf_list.txt

    # Add HLA_ALLELE annotation to each VCF before merging
    mkdir -p annotated_vcfs

    while IFS= read -r vcf_file; do
        # Extract allele name from filename
        # Assumes format: sample_ALLELE_caller.vcf.gz
        base=\$(basename "\$vcf_file" .vcf.gz)

        # Extract allele identifier from the filename
        # This regex extracts the allele part between sample name and caller
        allele=\$(echo "\$base" | sed -E 's/.*_([A-Z]_[0-9_]+)_.*/\\1/')

        if [ -n "\$allele" ] && [ "\$allele" != "\$base" ]; then
            # Add HLA_ALLELE to INFO field for each variant
            bcftools annotate \\
                --set-id '%CHROM\\_%POS\\_%REF\\_%ALT' \\
                --header-lines <(echo '##INFO=<ID=HLA_ALLELE,Number=1,Type=String,Description="HLA allele from which this variant was called">') \\
                "\$vcf_file" | \\
            bcftools annotate \\
                --annotations <(bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%CHROM_%POS_%REF_%ALT\\n' "\$vcf_file" | awk -v allele="\$allele" '{print \$1"\\t"\$2"\\t"\$5"\\t"allele}') \\
                --columns CHROM,POS,ID,INFO/HLA_ALLELE \\
                --output annotated_vcfs/\$(basename "\$vcf_file") \\
                --output-type z

            # Index the annotated VCF
            bcftools index -t annotated_vcfs/\$(basename "\$vcf_file")

            echo "Annotated \$vcf_file with allele: \$allele"
        else
            # If allele extraction fails, copy original
            echo "Warning: Could not extract allele from \$vcf_file, using original"
            cp "\$vcf_file" annotated_vcfs/
            cp "\${vcf_file}.tbi" annotated_vcfs/
        fi
    done < vcf_list.txt

    # Combine all annotated VCFs
    # Use concat if they're from different regions (alleles), merge if overlapping
    bcftools concat \\
        --allow-overlaps \\
        --remove-duplicates \\
        --output-type z \\
        --output ${prefix}.combined.vcf.gz \\
        ${args} \\
        annotated_vcfs/*.vcf.gz

    # Index the combined VCF
    bcftools index -t ${prefix}.combined.vcf.gz

    # Generate summary statistics
    echo "Combined VCF Statistics for ${meta.sample_id} (${caller}):" > ${prefix}.stats.txt
    echo "Total variants: \$(bcftools view -H ${prefix}.combined.vcf.gz | wc -l)" >> ${prefix}.stats.txt
    echo "" >> ${prefix}.stats.txt
    echo "Variants per allele:" >> ${prefix}.stats.txt
    bcftools query -f '%INFO/HLA_ALLELE\\n' ${prefix}.combined.vcf.gz | sort | uniq -c >> ${prefix}.stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
