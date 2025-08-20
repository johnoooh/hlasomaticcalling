process HLAHD {
    tag "${meta.id}"
    // label 'process_medium'

    // container = "cmopipeline/hlahd:1.4"
    container = "orgeraj/hlahd:1.7.1"
    scratch = true

    cpus = { 8 * task.attempt }
    memory = 5.GB

    // conda "bioconda::hlahd=1.7.0"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/hlahd:1.7.0--hdfd78af_0' :
    //     'biocontainers/hlahd:1.7.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${prefix}"), emit: results
    tuple val(meta), path("${prefix}_final.result.txt"), emit: hla_calls
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def fastq_files = fastq instanceof List ? fastq : [fastq]
    def fastq1 = fastq_files[0]
    def fastq2 = fastq_files.size() > 1 ? fastq_files[1] : ""
    """

    install_dir=/opt/hlahd/hlahd.1.7.1

    if [[ \$( ulimit -n ) -lt 1024 ]] ; then ulimit -n 1024 ;fi
        bash \$install_dir/bin/hlahd.sh -t ${task.cpus} -m 100 -f \$install_dir/freq_data \
        ${fastq1} \\
        ${fastq2} \\
        \$install_dir/HLA_gene.split.txt \
        \$install_dir/dictionary ${prefix} .

        cp */result/*_final.result.txt .


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hlahd: v1.7.1
    END_VERSIONS
    """
}