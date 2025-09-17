process SomaticCombineChannel {
  tag "CombineVCFS${meta.id}"
  container = "cmopipeline/bcftools-vt:1.2.3"

  // 3 intermidiate files (plus 3 index files) output for step by step filter check (2 filter steps involved here)
  // publishDir "${params.outDir}/somatic/${idTumor}__${idNormal}/combined_mutations/intermediate_files/", mode: params.publishDirMode, pattern: "*.union.annot.*"

  cpus = { 2 * task.attempt }
  memory = 5.GB

  input:
    tuple val(meta), path(mutectCombinedVcf), path(mutectCombinedVcfIndex), path(strelkaVcfSNV), path(strelkaVcfSNVIndex),path(strelkaVcfIndel), path(strelkaVcfIndelIndex)
    tuple val(meta2), path(fasta)
  output:
    tuple val(meta), path("*.union.vcf.gz"), emit: mutationMergedVcf
  
  when:
    task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"  

  """

  bcftools query -l ${mutectCombinedVcf} > samples.txt

  
  bcftools concat \
      --allow-overlaps \
      ${strelkaVcfSNV} ${strelkaVcfIndel} | \
  bcftools sort | \
  bcftools norm \
      --fasta-ref ${fasta} \
      --check-ref s \
      --output-type z \
      --output ${prefix}_strelka2.vcf.gz

  tabix --preset vcf ${prefix}_strelka2.vcf.gz


  merge_vcf.bash -n1 Mutect \\
  -n2 Strelka \\
  -f1 ${mutectCombinedVcf} \\
  -f2 ${prefix}_strelka2.vcf.gz \\
  -o ./tmp/ \\
  -p ${meta.id} \\
  -r1 samples.txt \\
  -r2 samples.txt 

  cp ./tmp/*.union.vcf.gz .

  """
}
