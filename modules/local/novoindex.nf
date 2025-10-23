nextflow.enable.dsl = 2

process NOVOINDEX {
    tag "$fasta"
    label 'process_high'

    conda "bioconda::novoalign=4.02.02"
    container "cmopipeline/lohhla:1.1.7"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.nix"), emit: index
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    novoindex \\
        ${args} \\
        ${prefix}.nix \\
        ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        novoalign: \$(novoindex 2>&1 | grep -oP 'V\\d+\\.\\d+\\.\\d+' | head -1 | sed 's/V//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.nix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        novoalign: \$(novoindex 2>&1 | grep -oP 'V\\d+\\.\\d+\\.\\d+' | head -1 | sed 's/V//')
    END_VERSIONS
    """
}