process BWA_MEM_CUSTOM {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bf/bf7890f8d4e38a7586581cb7fa13401b7af1582f21d94eef969df4cea852b6da/data' :
        'community.wave.seqera.io/library/bwa_htslib_samtools:56c9f8d5201889a4' }"

    input:
    tuple val(meta), path(reads), path(index), path(fasta)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Find the BWA index prefix

    INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`
    ls -l 
    echo \$INDEX
    bwa mem \\
        -k 23 \\
        -B 6  \\
        -O 8,8 \\
        -E 2,2 \\
        -L 10,10 \\
        -T 50 \\
        -t $task.cpus \\
        $args \\
        \$INDEX \\
        $reads \\
        | samtools $samtools_command $args2 ${reference} --threads $task.cpus -o ${prefix}.${extension} -

    samtools view -b -F 0x900 -f 0x2 -q 30 ${prefix}.${extension} -o aln.strict.filtered.bam
    samtools sort -@8 -o ${prefix}_aligned.bam aln.strict.filtered.bam 
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}