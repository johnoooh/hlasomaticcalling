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
include { GATK4_FILTERMUTECTCALLS } from '../modules/nf-core/gatk4/filtermutectcalls/main'
include { STRELKA_SOMATIC        } from '../modules/nf-core/strelka/somatic/main'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../modules/nf-core/gatk4/createsequencedictionary/main'
include { CREATE_HLA_REFERENCE } from '../modules/local/create_hla_reference'
include { BWA_MEM_CUSTOM } from '../modules/local/bwa_mem_custom'
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
    ch_reference = Channel.value(file(params.fasta, checkIfExists: true))
    ch_reference_fai = Channel.value(file(params.fastafai, checkIfExists: true))
    ch_hlahd_db = Channel.value(file(params.hlahd_db, checkIfExists: true))

    //
    // Parse input samplesheet - now expects normal_bam, normal_bai, tumor_bam, tumor_bai
    //
    // ch_samplesheet.view()
    ch_normal_bams = ch_samplesheet.map { row ->
        def meta = row[0]
        def normal_bam = row[1]
        def normal_bai = row[2]
        def normal_meta = meta + [sample_type: 'normal', id: "${meta.id}_normal"]
        [normal_meta, normal_bam, normal_bai]
    }
    // ch_samplesheet.view()
    ch_tumor_bams = ch_samplesheet.map { row ->
        def meta = row[0]
        def tumor_bam = row[3]
        def tumor_bai = row[4]
        def tumor_meta = meta + [sample_type: 'tumor', id: "${meta.id}_tumor"]
        [tumor_meta, tumor_bam, tumor_bai]
    }
    
    ch_all_bams = ch_normal_bams.mix(ch_tumor_bams)
    // ch_all_bams.view()

    //
    // Convert BAMs to FASTQ for processing
    //
    SAMTOOLS_FASTQ (
        ch_all_bams.map { meta, bam, bai -> [meta, bam] },
        false  // not interleaved
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions.first())

    // SAMTOOLS_FASTQ.out.fastq.view()
    //
    // MODULE: Run HLAHD on normal samples using FASTQ files
    //
    ch_normal_fastq = SAMTOOLS_FASTQ.out.fastq
        .filter { meta, fastq -> meta.sample_type == 'normal' }
    
    // ch_normal_fastq.view()
    HLAHD (
        ch_normal_fastq
    )

    ch_versions = ch_versions.mix(HLAHD.out.versions.first())

    
    //
    // Create HLA reference fastas and index them
    //
    ch_hla_calls = HLAHD.out.hla_calls
    // ch_hla_calls.view()

    CREATE_HLA_REFERENCE ( 
        ch_hla_calls, 
        ch_reference
    )
    // CREATE_HLA_REFERENCE.out.hla_reference.view()

    SAMTOOLS_FAIDX (
        CREATE_HLA_REFERENCE.out.hla_reference,
        CREATE_HLA_REFERENCE.out.hla_reference.map { meta, fasta -> [meta, []] },     // Channel 2: Empty fai input with correct meta
        true                                                                           // Channel 3: get_sizes parameter
    )
    // Create personalized HLA reference (this would need a custom process)
    // For now, we'll use the original reference and proceed with realignment
    
    BWA_INDEX (
        CREATE_HLA_REFERENCE.out.hla_reference,
    )

    ch_versions = ch_versions.mix(BWA_INDEX.out.versions.first())

    // BWA_INDEX.out.index.view { "BWA_INDEX output: $it" }
    // SAMTOOLS_FASTQ.out.fastq.view { "SAMTOOLS_FASTQ output: $it" }

    bwaindex    = BWA_INDEX.out.index

    GATK4_CREATESEQUENCEDICTIONARY (
        ch_reference.map { fasta -> [[id: 'reference'], fasta] }
    )

    ch_bwa_input = SAMTOOLS_FASTQ.out.fastq
    .combine(BWA_INDEX.out.index)
    .combine(CREATE_HLA_REFERENCE.out.hla_reference.map { fasta -> [[id: 'normal'], fasta] })

    ch_bwa_index = BWA_INDEX.out.index

    ch_fastq_with_patient = SAMTOOLS_FASTQ.out.fastq.map { meta, fastq ->
        def patient_id = meta.id.replace('_normal', '').replace('_tumor', '')
        [patient_id, meta, fastq]
    }

    ch_bwa_index_with_patient = BWA_INDEX.out.index.map { meta, index ->
        def patient_id = meta.id.replace('_normal', '')
        [patient_id, index]
    }
    
    ch_hla_ref_with_patient = CREATE_HLA_REFERENCE.out.hla_reference.map { meta, fasta ->
        def patient_id = meta.id.replace('_normal', '')
        [patient_id, fasta]
    }

    // Join all channels by patient ID to create complete alignment input
    ch_alignment_input = ch_fastq_with_patient
        .combine(ch_bwa_index_with_patient, by: 0)
        .combine(ch_hla_ref_with_patient, by: 0)
        .map { patient_id, sample_meta, fastq, index, fasta ->
            [sample_meta, fastq, index, fasta]
        }


    // ch_alignment_input.view()
    BWA_MEM (
        ch_alignment_input.map { sample_meta, fastq, index, fasta -> [sample_meta, fastq] },
        ch_alignment_input.map { sample_meta, fastq, index, fasta -> [['id': 'bwa_index'], index] },
        ch_alignment_input.map { sample_meta, fastq, index, fasta -> [['id': 'hla_reference'], fasta] },
        true
    )

    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    // Index realigned BAMs
    SAMTOOLS_INDEX (
        BWA_MEM.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    
    // Prepare tumor-normal pairs for somatic calling
    
    ch_realigned_bams = BWA_MEM.out.bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0])
    
    // ch_realigned_bams.view()

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

    // ch_tumor_normal_pairs.count().view { "Number of tumor-normal pairs: $it" }
    // ch_tumor_normal_pairs.view { "Tumor-normal pairs: $it" }

    //
    // MODULE: Run Mutect2 for somatic mutation calling
    //
    GATK4_MUTECT2 (
        ch_tumor_normal_pairs.map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
            [meta, [tumor_bam, normal_bam], [tumor_bai, normal_bai], []]
        },
        ch_reference.map { fasta -> [[id: 'reference'], fasta] },
        ch_reference_fai.map { fai -> [[id: 'reference'], fai] },
        GATK4_CREATESEQUENCEDICTIONARY.out.dict.map { meta, dict -> [[id: 'reference'], dict] },
        [],
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(GATK4_MUTECT2.out.versions.first())

    // Create the input channel for GATK4_FILTERMUTECTCALLS
    // This combines the VCF and stats outputs from MUTECT2
    // Create the input channel for GATK4_FILTERMUTECTCALLS
// This combines the VCF, TBI, and stats outputs from MUTECT2
    ch_filtermutect_in = GATK4_MUTECT2.out.vcf
    .join(GATK4_MUTECT2.out.tbi, by: 0)
    .join(GATK4_MUTECT2.out.stats, by: 0)
    .map { meta, vcf, tbi, stats ->
        [
            meta,           // meta map
            vcf,            // VCF file
            tbi,            // TBI index file
            stats,          // stats file
            [],             // orientationbias (empty - not using)
            [],             // segmentation (empty - not using)
            [],             // contamination table (empty - not using)
            0.0             // contamination estimate (0.0 - not using)
        ]
    }

GATK4_FILTERMUTECTCALLS(
    ch_filtermutect_in,
    ch_reference.map { fasta -> [[id: 'reference'], fasta] },
    ch_reference_fai.map { fai -> [[id: 'reference'], fai] },
    GATK4_CREATESEQUENCEDICTIONARY.out.dict.map { meta, dict -> [[id: 'reference'], dict] }
)
    //
    // MODULE: Run Strelka for somatic mutation calling
    //

    // Prepare personalized reference for Strelka
    ch_hla_ref_with_patient_for_strelka = CREATE_HLA_REFERENCE.out.hla_reference.map { meta, fasta ->
        def patient_id = meta.id.replace('_normal', '')
        [patient_id, fasta]
    }
    
    ch_hla_fai_with_patient_for_strelka = SAMTOOLS_FAIDX.out.fai.map { meta, fai ->
        def patient_id = meta.id.replace('_normal', '')
        [patient_id, fai]
    }
    // ch_hla_ref_with_patient_for_strelka.view()
    // ch_hla_fai_with_patient_for_strelka.view()


    ch_tumor_normal_keyed = ch_tumor_normal_pairs.map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
        [meta.id, meta, tumor_bam, tumor_bai, normal_bam, normal_bai]
    }

    // Join everything together
    ch_strelka_input = ch_tumor_normal_keyed
        .join(ch_hla_ref_with_patient_for_strelka)
        .join(ch_hla_fai_with_patient_for_strelka)
        .map { patient_id, meta, tumor_bam, tumor_bai, normal_bam, normal_bai, hla_fasta, hla_fai ->
            [meta, tumor_bam, tumor_bai, normal_bam, normal_bai, hla_fasta, hla_fai]
        }
    
    // ch_strelka_input.view()
    
    STRELKA_SOMATIC (
        ch_strelka_input.map { meta, normal_bam, normal_bai, tumor_bam, tumor_bai, hla_fasta, hla_fai -> [meta, normal_bam, normal_bai, tumor_bam, tumor_bai, [], [], [], []] },
        ch_strelka_input.map { meta, normal_bam, normal_bai, tumor_bam, tumor_bai, hla_fasta, hla_fai -> [meta, hla_fasta] },
        ch_strelka_input.map { meta, normal_bam, normal_bai, tumor_bam, tumor_bai, hla_fasta, hla_fai -> [meta, hla_fai] }
    
    )

    SomaticCombineMutect2Vcf.out.mutect2CombinedVcfOutput.combine(bamFiles, by: [0]).combine(SomaticRunStrelka2.out.strelka4Combine, by: [0,1,2]).set{ mutectStrelkaChannel }
    SomaticCombineChannel(mutectStrelkaChannel,
                          Channel.value([referenceMap.genomeFile, referenceMap.genomeIndex]),
                          Channel.value([referenceMap.repeatMasker, referenceMap.repeatMaskerIndex, referenceMap.mapabilityBlacklist, referenceMap.mapabilityBlacklistIndex]),
                          Channel.value([referenceMap.exomePoN, referenceMap.wgsPoN,referenceMap.exomePoNIndex, referenceMap.wgsPoNIndex,]),
                          Channel.value([referenceMap.gnomadWesVcf, referenceMap.gnomadWesVcfIndex,referenceMap.gnomadWgsVcf, referenceMap.gnomadWgsVcfIndex]))
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
    // mutect2_vcf    = GATK4_MUTECT2.out.vcf      // channel: Mutect2 VCF files
    // strelka_snvs   = STRELKA_SOMATIC.out.vcf_snvs    // channel: Strelka SNV VCF files
    // strelka_indels = STRELKA_SOMATIC.out.vcf_indels  // channel: Strelka indel VCF files

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
