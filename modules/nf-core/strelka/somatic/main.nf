process STRELKA_SOMATIC {
    tag "$meta.id"
    label 'process_medium'
    label 'error_retry'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/strelka:2.9.10--h9ee0642_1' :
        'biocontainers/strelka:2.9.10--h9ee0642_1' }"

    input:
    tuple val(meta), path(input_normal), path(input_index_normal), path(input_tumor), path(input_index_tumor),  path(manta_candidate_small_indels), path(manta_candidate_small_indels_tbi), path(target_bed), path(target_bed_index)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.somatic_indels.vcf.gz")    , emit: vcf_indels
    tuple val(meta), path("*.somatic_indels.vcf.gz.tbi"), emit: vcf_indels_tbi
    tuple val(meta), path("*.somatic_snvs.vcf.gz")      , emit: vcf_snvs
    tuple val(meta), path("*.somatic_snvs.vcf.gz.tbi")  , emit: vcf_snvs_tbi
    path "versions.yml"                                 , emit: versions
    tuple val(meta), path('*strelka2.vcf.gz'), path('*strelka2.vcf.gz.tbi'), emit: strelka4Combine

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def options_target_bed = target_bed ? "--callRegions ${target_bed}" : ""
    def options_manta = manta_candidate_small_indels ? "--indelCandidates ${manta_candidate_small_indels}" : ""
    """

    configureStrelkaSomaticWorkflow.py \\
        --tumor $input_tumor \\
        --normal $input_normal \\
        --referenceFasta $fasta \\
        ${options_target_bed} \\
        ${options_manta} \\
        $args \\
        --runDir strelka

    sed -i s/"isEmail = isLocalSmtp()"/"isEmail = False"/g strelka/runWorkflow.py

    python strelka/runWorkflow.py -m local -j $task.cpus
    mv strelka/results/variants/somatic.indels.vcf.gz     ${prefix}.somatic_indels.vcf.gz
    mv strelka/results/variants/somatic.indels.vcf.gz.tbi ${prefix}.somatic_indels.vcf.gz.tbi
    mv strelka/results/variants/somatic.snvs.vcf.gz       ${prefix}.somatic_snvs.vcf.gz
    mv strelka/results/variants/somatic.snvs.vcf.gz.tbi   ${prefix}.somatic_snvs.vcf.gz.tbi


    echo -e 'TUMOR ${prefix}_tumor\\nNORMAL ${prefix}_normal' > samples.txt

    bcftools concat \
        --allow-overlaps \
        ${prefix}.somatic_indels.vcf.gz ${prefix}.somatic_snvs.vcf.gz | \
    bcftools sort | \
    bcftools norm \
        --fasta-ref ${fasta} \
        --check-ref s \
        --output-type z \
        --output ${prefix}_strelka.vcf.gz

  tabix --preset vcf ${prefix}_strelka2.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strelka: \$( configureStrelkaSomaticWorkflow.py --version )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.somatic_indels.vcf.gz
    touch ${prefix}.somatic_indels.vcf.gz.tbi
    echo "" | gzip > ${prefix}.somatic_snvs.vcf.gz
    touch ${prefix}.somatic_snvs.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strelka: \$( configureStrelkaSomaticWorkflow.py --version )
    END_VERSIONS
    """
}
