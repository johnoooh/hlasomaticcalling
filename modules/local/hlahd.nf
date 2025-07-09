process HLAHD {
    tag "${meta.id}"
    label 'process_high'

    container = "cmopipeline/hlahd:1.4"
    scratch = true

    cpus = { 3 * task.attempt }
    memory = 10.GB

    // conda "bioconda::hlahd=1.7.0"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/hlahd:1.7.0--hdfd78af_0' :
    //     'biocontainers/hlahd:1.7.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(fastq)
    path hlahd_db
    path reference_fasta
    path reference_fai

    output:
    tuple val(meta), path("${prefix}"), emit: results
    tuple val(meta), path("${prefix}/result/${prefix}_final.result.txt"), emit: hla_calls
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
    # Convert BAM to FASTQ for HLAHD
    
    // # Run HLAHD
    // hlahd.sh \\
    //     -t ${task.cpus} \\
    //     -m 100 \\
    //     -c 0.95 \\
    //     ${args} \\
    //     ${prefix}_R1.fastq \\
    //     ${prefix}_R2.fastq \\
    //     ${hlahd_db}/HLA_gene.split.txt \\
    //     ${hlahd_db}/dictionary \\
    //     ${prefix} \\
    //     ${reference_fasta}

    install_dir=/hlahd.1.4.0

    if [[ \$( ulimit -n ) -lt 1024 ]] ; then ulimit -n 1024 ;fi
        bash \$install_dir/bin/hlahd.sh -t ${task.cpus} -m 100 -f \$install_dir/freq_data \
        ${fastq1} \\
        ${fastq2} \\
        \$install_dir/HLA_gene.split.txt \
        \$install_dir/dictionary ${prefix} .

        cp ${prefix}/result/${prefix}_final.result.txt .


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hlahd: \$(hlahd.sh 2>&1 | grep "HLA-HD version" | sed 's/.*version //g')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}