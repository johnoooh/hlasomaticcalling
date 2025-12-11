/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { HLAHD                  } from '../modules/local/hlahd'
include { BWA_INDEX              } from '../modules/nf-core/bwa/index/main'
include { BWA_MEM                } from '../modules/nf-core/bwa/mem/main'
include { NOVOALIGN              } from '../modules/local/novoalign'
include { NOVOINDEX              } from '../modules/local/novoindex'
include { SAMTOOLS_INDEX         } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FAIDX         } from '../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_FASTQ         } from '../modules/nf-core/samtools/fastq/main'
include { GATK4_MUTECT2          } from '../modules/nf-core/gatk4/mutect2/main'
include { GATK4_FILTERMUTECTCALLS } from '../modules/nf-core/gatk4/filtermutectcalls/main'
include { STRELKA_SOMATIC        } from '../modules/nf-core/strelka/somatic/main'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../modules/nf-core/gatk4/createsequencedictionary/main'
include { CREATE_HLA_REFERENCE } from '../modules/local/create_hla_reference'
include { BWA_MEM_CUSTOM } from '../modules/local/bwa_mem_custom'
include { PARSE_HLA_ALLELES } from '../modules/local/parse_hla_alleles'
include { EXTRACT_ALLELE_BAM } from '../modules/local/extract_allele_bam'
include { COMBINE_ALLELE_VCFS } from '../modules/local/combine_allele_vcfs'
include { EXTRACT_HLA_REGION } from '../modules/local/extract_hla_region'


include { SomaticCombineChannel } from '../modules/local/SomaticCombineChannel'
include { GENOMENEXUS_VCF2MAF } from '../modules/msk/genomenexus/vcf2maf/main'
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
    // MODULE: Extract HLA region from BAMs (POLYSOLVER approach)
    // This reduces data volume by extracting only HLA-relevant reads before FASTQ conversion
    //
    EXTRACT_HLA_REGION (
        ch_all_bams
    )
    ch_versions = ch_versions.mix(EXTRACT_HLA_REGION.out.versions.first())

    //
    // Convert HLA region BAMs to FASTQ for processing
    //
    SAMTOOLS_FASTQ (
        EXTRACT_HLA_REGION.out.bam.map { meta, bam, bai -> [meta, bam] },
        false
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
    NOVOINDEX(
        CREATE_HLA_REFERENCE.out.hla_reference,
    )

    // Create personalized HLA reference (this would need a custom process)
    // For now, we'll use the original reference and proceed with realignment
    
    // BWA_INDEX (
    //     CREATE_HLA_REFERENCE.out.hla_reference,
    // )

    // ch_versions = ch_versions.mix(BWA_INDEX.out.versions.first())

    // BWA_INDEX.out.index.view { "BWA_INDEX output: $it" }
    // SAMTOOLS_FASTQ.out.fastq.view { "SAMTOOLS_FASTQ output: $it" }

    novo_index    = NOVOINDEX.out.index

    GATK4_CREATESEQUENCEDICTIONARY (
        ch_reference.map { fasta -> [[id: 'reference'], fasta] }
    )


    ch_fastq_with_patient = SAMTOOLS_FASTQ.out.fastq.map { meta, fastq ->
        def patient_id = meta.id.replace('_normal', '').replace('_tumor', '')
        [patient_id, meta, fastq]
    }

    ch_novo_index_with_patient = NOVOINDEX.out.index.map { meta, index ->
        def patient_id = meta.id.replace('_normal', '')
        [patient_id, index]
    }
    
    ch_hla_ref_with_patient = CREATE_HLA_REFERENCE.out.hla_reference.map { meta, fasta ->
        def patient_id = meta.id.replace('_normal', '')
        [patient_id, fasta]
    }

    // Join all channels by patient ID to create complete alignment input
    ch_alignment_input = ch_fastq_with_patient
        .combine(ch_novo_index_with_patient, by: 0)
        .combine(ch_hla_ref_with_patient, by: 0)
        .map { patient_id, sample_meta, fastq, index, fasta ->
            [sample_meta, fastq, index, fasta]
        }


    // ch_alignment_input.view()
    NOVOALIGN (
        ch_alignment_input.map { sample_meta, fastq, index, fasta -> [sample_meta, fastq] },
        ch_alignment_input.map { sample_meta, fastq, index, fasta -> [['id': 'novo_index'], index] }
    )

    // ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    // BWA_MEM.out.bam.view()
    // Index realigned BAMs
    SAMTOOLS_INDEX (
        NOVOALIGN.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    //
    // MODULE: Parse HLA alleles for per-allele processing
    //
    PARSE_HLA_ALLELES (
        ch_hla_calls
    )
    ch_versions = ch_versions.mix(PARSE_HLA_ALLELES.out.versions.first())

    // Create a channel with individual alleles
    ch_alleles_per_sample = PARSE_HLA_ALLELES.out.alleles_list
        .flatMap { meta, alleles_file ->
            def alleles = alleles_file.readLines()
            alleles.collect { allele ->
                [meta.id.replace('_normal', ''), meta, allele.trim()]
            }
        }

    // Prepare tumor-normal pairs for somatic calling

    ch_realigned_bams = NOVOALIGN.out.bam
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

    // Prepare HLA reference for per-allele extraction
    ch_hla_ref_keyed = CREATE_HLA_REFERENCE.out.hla_reference
        .map { meta, fasta ->
            def patient_id = meta.id.replace('_normal', '')
            [patient_id, fasta]
        }

    // Combine tumor BAMs with alleles and HLA reference
    ch_tumor_for_extraction = ch_tumor_realigned
        .map { sample_id, meta, bam, bai -> [sample_id, meta, bam, bai] }
        .combine(ch_alleles_per_sample.map { sample_id, meta, allele -> [sample_id, allele] }, by: 0)
        .combine(ch_hla_ref_keyed, by: 0)
        .map { sample_id, meta, bam, bai, allele, hla_ref ->
            [meta, bam, bai, allele, hla_ref]
        }

    // Combine normal BAMs with alleles and HLA reference
    ch_normal_for_extraction = ch_normal_realigned
        .map { sample_id, meta, bam, bai -> [sample_id, meta, bam, bai] }
        .combine(ch_alleles_per_sample.map { sample_id, meta, allele -> [sample_id, allele] }, by: 0)
        .combine(ch_hla_ref_keyed, by: 0)
        .map { sample_id, meta, bam, bai, allele, hla_ref ->
            [meta, bam, bai, allele, hla_ref]
        }

    //
    // MODULE: Extract per-allele BAMs
    //
    EXTRACT_ALLELE_BAM (
        ch_tumor_for_extraction.mix(ch_normal_for_extraction)
    )
    ch_versions = ch_versions.mix(EXTRACT_ALLELE_BAM.out.versions.first())

    // Separate tumor and normal per-allele BAMs
    ch_tumor_allele_bams = EXTRACT_ALLELE_BAM.out.bam
        .filter { meta, bam, bai -> meta.sample_type == 'tumor' }
        .map { meta, bam, bai ->
            def sample_id = meta.id.replace('_tumor', '')
            [sample_id, meta.allele, meta, bam, bai]
        }

    ch_normal_allele_bams = EXTRACT_ALLELE_BAM.out.bam
        .filter { meta, bam, bai -> meta.sample_type == 'normal' }
        .map { meta, bam, bai ->
            def sample_id = meta.id.replace('_normal', '')
            [sample_id, meta.allele, meta, bam, bai]
        }

    // Join tumor and normal for each sample and allele
    ch_tumor_normal_pairs = ch_tumor_allele_bams
        .join(ch_normal_allele_bams, by: [0, 1])
        .map { sample_id, allele, tumor_meta, tumor_bam, tumor_bai, normal_meta, normal_bam, normal_bai ->
            def meta = [
                id: "${sample_id}_${tumor_meta.allele_safe}_somatic",
                sample_id: sample_id,
                allele: allele,
                allele_safe: tumor_meta.allele_safe,
                tumor_id: tumor_meta.id,
                normal_id: normal_meta.id
            ]
            [meta, tumor_bam, tumor_bai, normal_bam, normal_bai]
        }

    ch_tumor_normal_pairs.count().view { "Number of tumor-normal-allele pairs: $it" }
    // ch_tumor_normal_pairs.view { "Tumor-normal-allele pairs: $it" }

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
    // Create a channel with allele counts per sample for incremental grouping
    //
    ch_allele_counts = PARSE_HLA_ALLELES.out.alleles_list
        .map { meta, alleles_file ->
            def sample_id = meta.id.replace('_normal', '')
            def alleles = alleles_file.readLines()
            [sample_id, alleles.size()]
        }

    ch_allele_counts.view { "DEBUG ch_allele_counts: $it" }

    //
    // Combine per-allele VCFs into per-sample VCFs (Mutect2)
    //
    ch_mutect_before_combine = GATK4_FILTERMUTECTCALLS.out.vcf
        .join(GATK4_FILTERMUTECTCALLS.out.tbi, by: 0)
        .map { meta, vcf, tbi ->
            def sample_id = meta.sample_id
            [sample_id, meta, vcf, tbi]
        }

    ch_mutect_before_combine.view { "DEBUG ch_mutect_before_combine: sample_id=${it[0]}" }

    ch_mutect_per_sample = ch_mutect_before_combine
        .combine(ch_allele_counts, by: 0)
        .map { sample_id, meta, vcf, tbi, count ->
            // Include count in grouping key for independent completion
            [sample_id, 'mutect2', count, meta, vcf, tbi]
        }
        .groupTuple(by: [0, 1, 2])  // Group by sample_id, caller, AND expected count
        .map { sample_id, caller, count, metas, vcfs, tbis ->
            // Debug: show what we got
            println "DEBUG grouped mutect: sample=${sample_id}, count=${count}, vcfs.size=${vcfs.size()}"
            [sample_id, caller, count, metas, vcfs, tbis]
        }
        .filter { sample_id, caller, count, metas, vcfs, tbis ->
            // Only emit when we have all alleles for this sample
            def pass = (vcfs.size() == count)
            if (!pass) {
                println "DEBUG filter blocked: sample=${sample_id}, expected=${count}, got=${vcfs.size()}"
            }
            return pass
        }
        .map { sample_id, caller, count, metas, vcfs, tbis ->
            def meta = [sample_id: sample_id, id: sample_id]
            [meta, caller, vcfs, tbis]
        }

    //
    // MODULE: Run Strelka for somatic mutation calling
    //
    // CREATE_HLA_REFERENCE.out.hla_reference.view()
    // Prepare personalized reference for Strelka (per-allele)
    // The HLA reference contains all alleles, but BAMs are already filtered per-allele
    ch_hla_ref_for_strelka = CREATE_HLA_REFERENCE.out.hla_reference
        .map { meta, fasta ->
            def patient_id = meta.id.replace('_normal', '')
            [patient_id, fasta]
        }

    ch_hla_fai_for_strelka = SAMTOOLS_FAIDX.out.fai
        .map { meta, fai ->
            def patient_id = meta.id.replace('_normal', '')
            [patient_id, fai]
        }

    // Key tumor-normal pairs by sample_id for joining with references
    ch_tumor_normal_keyed = ch_tumor_normal_pairs.map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai ->
        [meta.sample_id, meta, tumor_bam, tumor_bai, normal_bam, normal_bai]
    }

    // Join everything together
    ch_strelka_input = ch_tumor_normal_keyed
        .combine(ch_hla_ref_for_strelka, by: 0)
        .combine(ch_hla_fai_for_strelka, by: 0)
        .map { patient_id, meta, tumor_bam, tumor_bai, normal_bam, normal_bai, hla_fasta, hla_fai ->
            [meta, normal_bam, normal_bai, tumor_bam, tumor_bai, hla_fasta, hla_fai]
        }
    
    // ch_strelka_input.view()
    
    STRELKA_SOMATIC (
        ch_strelka_input.map { meta, normal_bam, normal_bai, tumor_bam, tumor_bai, hla_fasta, hla_fai -> [meta, normal_bam, normal_bai, tumor_bam, tumor_bai, [], [], [], []] },
        ch_strelka_input.map { meta, normal_bam, normal_bai, tumor_bam, tumor_bai, hla_fasta, hla_fai -> [meta, hla_fasta] },
        ch_strelka_input.map { meta, normal_bam, normal_bai, tumor_bam, tumor_bai, hla_fasta, hla_fai -> [meta, hla_fai] }

    )

    ch_versions = ch_versions.mix(STRELKA_SOMATIC.out.versions.first())

    //
    // Combine per-allele VCFs into per-sample VCFs (Strelka SNVs)
    //
    ch_strelka_snvs_per_sample = STRELKA_SOMATIC.out.vcf_snvs
        .join(STRELKA_SOMATIC.out.vcf_snvs_tbi, by: 0)
        .map { meta, vcf, tbi ->
            def sample_id = meta.sample_id
            [sample_id, meta, vcf, tbi]
        }
        .combine(ch_allele_counts, by: 0)
        .map { sample_id, meta, vcf, tbi, count ->
            [sample_id, 'strelka_snvs', count, meta, vcf, tbi]
        }
        .groupTuple(by: [0, 1, 2])
        .filter { sample_id, caller, count, metas, vcfs, tbis ->
            vcfs.size() == count
        }
        .map { sample_id, caller, count, metas, vcfs, tbis ->
            def meta = [sample_id: sample_id, id: sample_id]
            [meta, caller, vcfs, tbis]
        }

    //
    // Combine per-allele VCFs into per-sample VCFs (Strelka Indels)
    //
    ch_strelka_indels_per_sample = STRELKA_SOMATIC.out.vcf_indels
        .join(STRELKA_SOMATIC.out.vcf_indels_tbi, by: 0)
        .map { meta, vcf, tbi ->
            def sample_id = meta.sample_id
            [sample_id, meta, vcf, tbi]
        }
        .combine(ch_allele_counts, by: 0)
        .map { sample_id, meta, vcf, tbi, count ->
            [sample_id, 'strelka_indels', count, meta, vcf, tbi]
        }
        .groupTuple(by: [0, 1, 2])
        .filter { sample_id, caller, count, metas, vcfs, tbis ->
            vcfs.size() == count
        }
        .map { sample_id, caller, count, metas, vcfs, tbis ->
            def meta = [sample_id: sample_id, id: sample_id]
            [meta, caller, vcfs, tbis]
        }

    //
    // MODULE: Combine per-allele VCFs into per-sample VCFs
    //
    COMBINE_ALLELE_VCFS (
        ch_mutect_per_sample.mix(
            ch_strelka_snvs_per_sample,
            ch_strelka_indels_per_sample
        )
    )
    ch_versions = ch_versions.mix(COMBINE_ALLELE_VCFS.out.versions.first())

    //
    // Prepare combined VCFs for SomaticCombineChannel
    // Join Mutect2, Strelka SNVs, and Strelka Indels by sample_id
    //
    ch_combined_mutect = COMBINE_ALLELE_VCFS.out.vcf
        .filter { meta, vcf, tbi -> vcf.name.contains('mutect2') }
        .map { meta, vcf, tbi -> [meta.sample_id, meta, vcf, tbi] }

    ch_combined_strelka_snvs = COMBINE_ALLELE_VCFS.out.vcf
        .filter { meta, vcf, tbi -> vcf.name.contains('strelka_snvs') }
        .map { meta, vcf, tbi -> [meta.sample_id, meta, vcf, tbi] }

    ch_combined_strelka_indels = COMBINE_ALLELE_VCFS.out.vcf
        .filter { meta, vcf, tbi -> vcf.name.contains('strelka_indels') }
        .map { meta, vcf, tbi -> [meta.sample_id, meta, vcf, tbi] }

    // Join all three by sample_id
    ch_for_somatic_combine = ch_combined_mutect
        .join(ch_combined_strelka_snvs, by: 0)
        .join(ch_combined_strelka_indels, by: 0)
        .map { sample_id, mutect_meta, mutect_vcf, mutect_tbi,
               strelka_snv_meta, strelka_snv_vcf, strelka_snv_tbi,
               strelka_indel_meta, strelka_indel_vcf, strelka_indel_tbi ->
            // Use mutect_meta as the base meta
            [mutect_meta, mutect_vcf, mutect_tbi,
             strelka_snv_vcf, strelka_snv_tbi,
             strelka_indel_vcf, strelka_indel_tbi]
        }

    SomaticCombineChannel(
        ch_for_somatic_combine,
        ch_reference.map { fasta -> [[id: 'reference'], fasta] }
    )


    // GENOMENEXUS_VCF2MAF(SomaticCombineChannel.out.mutationMergedVcf)



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
    hla_region_stats = EXTRACT_HLA_REGION.out.stats // channel: HLA region extraction statistics
    hla_calls      = HLAHD.out.hla_calls        // channel: HLA typing results
    alleles_list   = PARSE_HLA_ALLELES.out.alleles_list // channel: Per-sample allele list
    mutect2_vcf    = GATK4_MUTECT2.out.vcf      // channel: Per-allele Mutect2 VCF files
    strelka_snvs   = STRELKA_SOMATIC.out.vcf_snvs    // channel: Per-allele Strelka SNV VCF files
    strelka_indels = STRELKA_SOMATIC.out.vcf_indels  // channel: Per-allele Strelka indel VCF files
    allele_bams    = EXTRACT_ALLELE_BAM.out.bam     // channel: Per-allele BAM files
    combined_vcfs  = COMBINE_ALLELE_VCFS.out.vcf    // channel: Per-sample combined VCF files

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
