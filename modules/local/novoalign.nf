#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process NOVOALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::novoalign=4.02.02 bioconda::samtools=1.17"
    container "cmopipeline/lohhla:1.1.7"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*.sorted.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def samtools_args = task.ext.samtools_args ?: '-b -h'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def forward_reads = reads[0]  // First file is the forward reads
    def reverse_reads = reads[1]  // Second file is the reverse reads
   
    def paired_end = meta.single_end ? '' : '-i PE 250,30'
    """


    # Uncompress the input files to the temporary directory
    gunzip -c ${forward_reads} > forward_reads.fastq
    gunzip -c ${reverse_reads} > reverse_reads.fastq


    # Run Novoalign using the uncompressed file
    novoalign ${args} -d ${index} -f forward_reads.fastq reverse_reads.fastq ${paired_end} -o SAM ${args2} | samtools view --threads ${task.cpus} ${samtools_args} -o ${prefix}.bam -


    # Clean up the temporary directory
    samtools sort -o ${prefix}.sorted.bam ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        novoalign: \$(novoalign 2>&1 | grep -oP 'V\\d+\\.\\d+\\.\\d+' | head -1 | sed 's/V//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        novoalign: \$(novoalign 2>&1 | grep -oP 'V\\d+\\.\\d+\\.\\d+' | head -1 | sed 's/V//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}