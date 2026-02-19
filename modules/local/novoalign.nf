#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process NOVOALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::novoalign=4.02.02 bioconda::samtools=1.17"
    container "docker.io/cmopipeline/lohhla:1.1.7"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*.sorted.bam"), emit: bam
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def samtools_args = task.ext.samtools_args ?: '-b -h'
    def prefix = task.ext.prefix ?: "${meta.id}"
    def forward_reads = reads[0]
    def reverse_reads = reads[1]

    // Polysolver uses -o FullNW by default (no soft clipping)
    // Set soft_clip parameter via task.ext.soft_clip if you want soft clipping mode
    def soft_clip = task.ext.soft_clip ?: false
    def alignment_mode = soft_clip ? '-g 20 -x 3' : '-o FullNW'

    // Prepare read group string for samtools addreplacerg
    def read_group = "ID:${meta.id}\\tSM:${meta.id}\\tPL:ILLUMINA\\tLB:${meta.id}\\tPU:${meta.id}"

    """
    # Uncompress the input files
    gunzip -c ${forward_reads} > forward_reads.fastq
    gunzip -c ${reverse_reads} > reverse_reads.fastq

    # Run Novoalign with Polysolver parameters
    # -R 0: Report threshold 0 (all alignments)
    # -r all: Report ALL alignments meeting threshold (critical for HLA multi-mapping)
    # -o SAM: SAM output format
    # -o FullNW OR -g 20 -x 3: Alignment mode (no soft-clip vs soft-clip)
    # grep -P '\\thla': Keep only HLA-aligned reads
    # samtools addreplacerg: Add read group to BAM
    novoalign -d ${index} \\
        -f forward_reads.fastq reverse_reads.fastq \\
        -F STDFQ \\
        ${args} \\
        ${alignment_mode} \\
        -o SAM | \\
        samtools view --threads ${task.cpus} ${samtools_args} - | \\
        samtools addreplacerg -r '@RG\\t${read_group}' -o ${prefix}.bam -

    # Sort the BAM file
    samtools sort --threads ${task.cpus} -o ${prefix}.sorted.bam ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        novoalign: \$(novoalign 2>&1 | grep -oP 'V\\d+\\.\\d+\\.\\d+' | head -1 | sed 's/V//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        novoalign: \$(novoalign 2>&1 | grep -oP 'V\\d+\\.\\d+\\.\\d+' | head -1 | sed 's/V//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}