process HAPSTER {
    tag "${meta.id}"
    label 'process_high'

    container = "orgeraj/hapster:1.0.1"
    scratch = true

    cpus = { 8 * task.attempt }
    memory = { 48.GB * task.attempt }
    time = { 24.h * task.attempt }

    input:
    tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
    tuple val(meta2), path(hla_calls)
    path(hapster_refs)
    path(extraction_regions)

    output:
    tuple val(meta), path("${prefix}/somatic/*.somatic.filtered.vcf.gz"),     emit: somatic_vcf
    tuple val(meta), path("${prefix}/somatic/*.somatic.filtered.vcf.gz.tbi"), emit: somatic_vcf_tbi
    tuple val(meta), path("${prefix}/somatic/*.kmer_filtered.vcf.gz"),        emit: somatic_vcf_kmer_filtered
    tuple val(meta), path("${prefix}/somatic/*.kmer_filtered.vcf.gz.tbi"),    emit: somatic_vcf_kmer_filtered_tbi
    tuple val(meta), path("${prefix}/somatic/*.annotated.txt"),               emit: annotated_variants, optional: true
    tuple val(meta), path("${prefix}/germline/*.germline.vcf.gz"),            emit: germline_vcf, optional: true
    tuple val(meta), path("${prefix}/germline/*.germline.vcf.gz.tbi"),        emit: germline_vcf_tbi, optional: true
    tuple val(meta), path("${prefix}/haplotypes"),                            emit: haplotype_dir
    path "versions.yml",                                                      emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def protocol = task.ext.protocol ?: 'wgs'
    def genes = task.ext.genes ?: 'A,B,C'
    def sample_id = meta.id
    def tumor_id = meta.tumor_id ?: "${meta.id}_tumor"
    def normal_id = meta.normal_id ?: "${meta.id}_normal"
    """
    # Step 1: Convert HLAHD calls to Hapster haplotype format
    mkdir -p ${prefix}/haplotypes
    convert_hlahd_to_hapster.py \\
        --hlahd_results ${hla_calls} \\
        --output_dir ${prefix}/haplotypes \\
        --genes ${genes}

    # Step 2: Extract HLA reads (normal)
    hapster extract_dna \\
        --bam ${normal_bam} \\
        --sample-name ${normal_id} \\
        --output-dir ${prefix}/extract_normal \\
        --ref-dir ${hapster_refs} \\
        --regions ${extraction_regions} \\
        --threads ${task.cpus} \\
        ${args}

    # Step 3: Extract HLA reads (tumor)
    hapster extract_dna \\
        --bam ${tumor_bam} \\
        --sample-name ${tumor_id} \\
        --output-dir ${prefix}/extract_tumor \\
        --ref-dir ${hapster_refs} \\
        --regions ${extraction_regions} \\
        --threads ${task.cpus} \\
        ${args}

    # Step 4: Call germline mutations against personalized reference
    hapster germline_mutations \\
        --extracted-dir ${prefix}/extract_normal \\
        --haplotype-file ${prefix}/haplotypes/haplotypes.csv \\
        --output-dir ${prefix}/germline \\
        --ref-dir ${hapster_refs} \\
        --threads ${task.cpus} \\
        --protocol ${protocol}

    # Step 5: Realign normal reads to germline-imputed reference
    hapster dna_realign \\
        --extracted-dir ${prefix}/extract_normal \\
        --germline-dir ${prefix}/germline \\
        --output-dir ${prefix}/realign_normal \\
        --ref-dir ${hapster_refs} \\
        --threads ${task.cpus}

    # Step 6: Realign tumor reads to germline-imputed reference
    hapster dna_realign \\
        --extracted-dir ${prefix}/extract_tumor \\
        --germline-dir ${prefix}/germline \\
        --output-dir ${prefix}/realign_tumor \\
        --ref-dir ${hapster_refs} \\
        --threads ${task.cpus}

    # Step 7: Somatic mutation calling (Mutect2 + kmer filtering)
    hapster somatic_mutations \\
        --normal-dir ${prefix}/realign_normal \\
        --tumor-dir ${prefix}/realign_tumor \\
        --germline-dir ${prefix}/germline \\
        --output-dir ${prefix}/somatic \\
        --ref-dir ${hapster_refs} \\
        --threads ${task.cpus} \\
        --normal-name ${normal_id} \\
        --tumor-name ${tumor_id}

    # Step 8: Bgzip and tabix any uncompressed VCFs
    find ${prefix}/somatic -name "*.vcf" -not -name "*.vcf.gz" | while read vcf; do
        bgzip -c "\${vcf}" > "\${vcf}.gz"
        tabix -p vcf "\${vcf}.gz"
    done

    find ${prefix}/germline -name "*.vcf" -not -name "*.vcf.gz" | while read vcf; do
        bgzip -c "\${vcf}" > "\${vcf}.gz"
        tabix -p vcf "\${vcf}.gz"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hapster: 1.0.1
        convert_hlahd_to_hapster: 1.0.0
    END_VERSIONS
    """
}
