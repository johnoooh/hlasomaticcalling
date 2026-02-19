#!/usr/bin/env nextflow

/*
 * Standalone test: Verify CREATE_HLA_REFERENCE works with the new
 * IMGT-based decoy logic (no genome_fasta dependency).
 *
 * Usage:
 *   nextflow run test_create_hla_ref.nf -profile docker --outdir test_output
 */

nextflow.enable.dsl = 2

include { CREATE_HLA_REFERENCE } from './modules/local/create_hla_reference'

params.outdir = 'test_output'
params.fasta = "${projectDir}/assets/imgt_hla_gen.fasta.gz"
params.fastafai = "${projectDir}/assets/imgt_hla_gen.fasta.gz.fai"

workflow {
    // Mock HLA typing calls (simulates HLAHD output)
    ch_hla_calls = Channel.of(
        [ [id: 'test_normal'], file("${projectDir}/test_data/mock_hla_calls.txt") ]
    )

    ch_reference = Channel.value(file(params.fasta, checkIfExists: true))

    CREATE_HLA_REFERENCE(
        ch_hla_calls,
        ch_reference
    )

    // Verify the output
    CREATE_HLA_REFERENCE.out.hla_reference
        .map { meta, fasta ->
            println "=== CREATE_HLA_REFERENCE output ==="
            println "Meta: ${meta}"
            println "FASTA: ${fasta}"
            println "File size: ${fasta.size()} bytes"
            [meta, fasta]
        }
        .set { ch_result }

    // Run validation
    VALIDATE_REFERENCE(ch_result)
}

process VALIDATE_REFERENCE {
    tag "${meta.id}"
    publishDir "${params.outdir}", mode: 'copy'

    conda "bioconda::samtools=1.19.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0' :
        'biocontainers/samtools:1.19.2--h50ea8bc_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    path "validation_report.txt"

    script:
    """
    echo "=== HLA Reference Validation ===" > validation_report.txt
    echo "" >> validation_report.txt

    echo "--- All contigs ---" >> validation_report.txt
    grep "^>" ${fasta} >> validation_report.txt
    echo "" >> validation_report.txt

    echo "--- Patient alleles (A/B/C, no _decoy suffix) ---" >> validation_report.txt
    grep "^>" ${fasta} | grep -v "_decoy" >> validation_report.txt || echo "(none)" >> validation_report.txt
    echo "" >> validation_report.txt

    echo "--- Decoy contigs ---" >> validation_report.txt
    grep "^>" ${fasta} | grep "_decoy" >> validation_report.txt || echo "(none)" >> validation_report.txt
    echo "" >> validation_report.txt

    echo "--- Decoy gene summary ---" >> validation_report.txt
    grep "^>" ${fasta} | grep "_decoy" | sed 's/>\\([^_]*\\)_.*/\\1/' | sort -u >> validation_report.txt
    echo "" >> validation_report.txt

    TOTAL=\$(grep -c "^>" ${fasta})
    PATIENT=\$(grep "^>" ${fasta} | grep -cv "_decoy" || echo 0)
    DECOY=\$(grep "^>" ${fasta} | grep -c "_decoy" || echo 0)
    echo "--- Counts ---" >> validation_report.txt
    echo "Total contigs: \$TOTAL" >> validation_report.txt
    echo "Patient alleles: \$PATIENT" >> validation_report.txt
    echo "Decoy contigs: \$DECOY" >> validation_report.txt

    # Verify no A/B/C decoys (those are patient-typed)
    echo "" >> validation_report.txt
    echo "--- Checks ---" >> validation_report.txt
    if grep "^>" ${fasta} | grep "_decoy" | grep -qE "^>A_|^>B_|^>C_"; then
        echo "FAIL: Found A/B/C decoy contigs (should be excluded)" >> validation_report.txt
    else
        echo "PASS: No A/B/C decoy contigs found" >> validation_report.txt
    fi

    # Verify DRB1 decoy exists (Class II coverage)
    if grep "^>" ${fasta} | grep -q "DRB1.*_decoy"; then
        echo "PASS: DRB1 decoy present (Class II coverage)" >> validation_report.txt
    else
        echo "FAIL: No DRB1 decoy found" >> validation_report.txt
    fi

    # Verify E decoy exists (non-classical Class I)
    if grep "^>" ${fasta} | grep -q "E_.*_decoy"; then
        echo "PASS: HLA-E decoy present (non-classical Class I)" >> validation_report.txt
    else
        echo "FAIL: No HLA-E decoy found" >> validation_report.txt
    fi

    # Verify Y decoy exists (was previously hardcoded CDS)
    if grep "^>" ${fasta} | grep -q "Y_.*_decoy"; then
        echo "PASS: HLA-Y decoy present (previously hardcoded)" >> validation_report.txt
    else
        echo "FAIL: No HLA-Y decoy found" >> validation_report.txt
    fi

    # Verify MICA decoy exists (new coverage)
    if grep "^>" ${fasta} | grep -q "MICA.*_decoy"; then
        echo "PASS: MICA decoy present (new coverage)" >> validation_report.txt
    else
        echo "FAIL: No MICA decoy found" >> validation_report.txt
    fi

    cat validation_report.txt
    """
}
