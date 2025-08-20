process SomaticCombineChannel {
  tag "${idTumor + "__" + idNormal}"

  // 3 intermidiate files (plus 3 index files) output for step by step filter check (2 filter steps involved here)
  // publishDir "${params.outDir}/somatic/${idTumor}__${idNormal}/combined_mutations/intermediate_files/", mode: params.publishDirMode, pattern: "*.union.annot.*"

  input:
    tuple val(meta), path(mutectCombinedVcf), path(mutectCombinedVcfIndex), path(strelkaVcf), path(strelkaVcfIndex)

  output:
    tuple val(idTumor), val(idNormal), val(target), path("${outputPrefix}.pass.vcf"), emit: mutationMergedVcf
    path("${idTumor}__${idNormal}.union.annot.vcf.gz")
    path("${idTumor}__${idNormal}.union.annot.vcf.gz.tbi")
    path("${idTumor}__${idNormal}.union.annot.filter.vcf.gz")
    path("${idTumor}__${idNormal}.union.annot.filter.vcf.gz.tbi")
    path("${idTumor}__${idNormal}.union.annot.filter.pass.vcf.gz")
    path("${idTumor}__${idNormal}.union.annot.filter.pass.vcf.gz.tbi")
  
  script:
  outputPrefix = "${idTumor}__${idNormal}"
  isecDir = "${idTumor}.isec"
  pon = wgsPoN
  gnomad = gnomadWgsVcf
  if (target == "wgs") {
    pon = wgsPoN
    gnomad = gnomadWgsVcf
  }
  else {
    pon = exomePoN
    gnomad = gnomadWesVcf
  }
  """

  echo -e 'TUMOR ${idTumor}\\nNORMAL ${idNormal}' > samples.txt
  
  bcftools concat \
    --allow-overlaps \
    Strelka_${outputPrefix}_somatic_indels.vcf.gz Strelka_${outputPrefix}_somatic_snvs.vcf.gz | \
  bcftools reheader \
    --samples samples.txt | \
  bcftools sort | \
  bcftools norm \
    --fasta-ref ${genomeFile} \
    --check-ref s \
    --output-type z \
    --output ${outfile}

  tabix --preset vcf ${outfile}


  merge_vcf.bash -n1 Mutect \\
  -n2 Strelka \\
  -f1 ${mutectCombinedVcf} \\
  -f2 ${strelkaVcf} \\
  -o outtmp
  -p ${meta_id}
  -r1 
  """
}
