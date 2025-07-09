/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { HLAHD                  } from '../modules/local/hlahd'
include { BWA_INDEX              } from '../modules/nf-core/bwa/index/main'
include { BWA_MEM                } from '../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_INDEX         } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FAIDX         } from '../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FASTQ         } from '../modules/nf-core/samtools/fastq/main'
include { GATK4_MUTECT2          } from '../modules/nf-core/gatk4/mutect2/main'
include { STRELKA_SOMATIC        } from '../modules/nf-core/strelka/somatic/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_hlasomatic_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow HLASOMATIC {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // Prepare reference files
    //
    ch_reference = Channel.fromPath(params.fasta, checkIfExists: true)
    ch_hlahd_db = Channel.fromPath(params.hlahd_db, checkIfExists: true)
    
    // Index reference fasta
    SAMTOOLS_FAIDX (
        ch_reference.map { fasta -> [[id: 'reference'], fasta] }
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())

    //
    // Parse input samplesheet - now expects normal_bam, normal_bai, tumor_bam, tumor_bai
    //
    ch_normal_bams = ch_samplesheet.map { meta, normal_bam, normal_bai, tumor_bam, tumor_bai ->
        def normal_meta = meta + [sample_type: 'normal', id: "${meta.id}_normal"]
        [normal_meta, normal_bam, normal_bai]
    }
    
    ch_tumor_bams = ch_samplesheet.map { meta, normal_bam, normal_bai, tumor_bam, tumor_bai ->
        def tumor_meta = meta + [sample_type: 'tumor', id: "${meta.id}_tumor"]
        [tumor_meta, tumor_bam, tumor_bai]
    }
    
    ch_all_bams = ch_normal_bams.mix(ch_tumor_bams)

    //
    // Convert BAMs to FASTQ for processing
    //
    SAMTOOLS_FASTQ (
        ch_all_bams.map { meta, bam, bai -> [meta, bam] },
        false  // not interleaved
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions.first())

    //
    // MODULE: Run HLAHD on normal samples using FASTQ files
    //
    ch_normal_fastq = SAMTOOLS_FASTQ.out.fastq
        .filter { meta, fastq -> meta.sample_type == 'normal' }
    
    HLAHD (
        ch_normal_fastq,
        ch_hlahd_db,
        ch_reference,
        SAMTOOLS_FAIDX.out.fai.map { meta, fai -> fai }
    )
    ch_versions = ch_versions.mix(HLAHD.out.versions.first())

    //
    // Create HLA reference fastas and index them
    //
    ch_hla_calls = HLAHD.out.hla_calls
    
    // Create personalized HLA reference (this would need a custom process)
    // For now, we'll use the original reference and proceed with realignment
    
    BWA_INDEX (
        ch_reference.map { fasta -> [[id: 'reference'], fasta] }
    )
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions.first())

    //
    // Realign both tumor and normal samples using converted FASTQ files
    //
    BWA_MEM (
        SAMTOOLS_FASTQ.out.fastq,
        BWA_INDEX.out.index.map { meta, index -> index },
        ch_reference.map { fasta -> [[id: 'reference'], fasta] },
        true
    )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    // Index realigned BAMs
    SAMTOOLS_INDEX (
        BWA_MEM.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    //
    // Prepare tumor-normal pairs for somatic calling
    //
    ch_realigned_bams = BWA_MEM.out.bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0])
    
    // Get tumor realigned BAMs
    ch_tumor_realigned = ch_realigned_bams
        .filter { meta, bam, bai -> meta.sample_type == 'tumor' }
        .map { meta, bam, bai -> [meta.id.replace('_tumor', ''), meta, bam, bai] }
    
    // Get normal realigned BAMs  
    ch_normal_realigned = ch_realigned_bams
        .filter { meta, bam, bai -> meta.sample_type == 'normal' }
        .map { meta, bam, bai -> [meta.id.replace('_normal', ''), meta, bam, bai] }
    
    // Join tumor and normal for each sample
    ch_tumor_normal_pairs = ch_tumor_realigned
        .join(ch_normal_realigned)
        .map { sample_id, tumor_meta, tumor_bam, tumor_bai, normal_meta, normal_bam, normal_bai ->
            def meta = [
                id: "${sample_id}_somatic",
                sample_id: sample_id,
                tumor_id: tumor_meta.id,
                normal_id: normal_meta.id
            ]
            [meta, tumor_bam, tumor_bai, normal_bam, normal_bai]
        }

    //
    // MODULE: Run Mutect2 for somatic mutation calling
    //
    GATK4_MUTECT2 (
        ch_tumor_normal_pairs.map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
            [meta, [tumor_bam, normal_bam], [tumor_bai, normal_bai], []]
        },
        ch_reference.map { fasta -> [[id: 'reference'], fasta] },
        SAMTOOLS_FAIDX.out.fai.map { meta, fai -> [[id: 'reference'], fai] },
        [],
        [],
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(GATK4_MUTECT2.out.versions.first())

    //
    // MODULE: Run Strelka for somatic mutation calling
    //
    STRELKA_SOMATIC (
        ch_tumor_normal_pairs.map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
            [meta, normal_bam, normal_bai, tumor_bam, tumor_bai, [], [], [], []]
        },
        ch_reference,
        SAMTOOLS_FAIDX.out.fai.map { meta, fai -> fai }
    )
    ch_versions = ch_versions.mix(STRELKA_SOMATIC.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'hlasomatic_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
    hla_calls      = HLAHD.out.hla_calls        // channel: HLA typing results
    mutect2_vcf    = GATK4_MUTECT2.out.vcf      // channel: Mutect2 VCF files
    strelka_snvs   = STRELKA_SOMATIC.out.vcf_snvs    // channel: Strelka SNV VCF files
    strelka_indels = STRELKA_SOMATIC.out.vcf_indels  // channel: Strelka indel VCF files

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
