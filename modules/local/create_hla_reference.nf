process CREATE_HLA_REFERENCE {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::biopython=1.81 bioconda::pysam=0.22.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-3a59640f3fe1ed11819984087d31d68600200c3f:185a25ca79923df85b58f42deb48f5ac4481e91f-0' :
        'biocontainers/mulled-v2-3a59640f3fe1ed11819984087d31d68600200c3f:185a25ca79923df85b58f42deb48f5ac4481e91f-0' }"
    input:
    tuple val(meta), path(hla_calls)
    path reference_fasta
    path hla_sequences_db

    output:
    tuple val(meta), path("${prefix}_hla_reference.fasta"), emit: hla_reference
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3
    
    import pysam 
    # For now, just copy the original reference
    
    geneset = set(["A", B", "C" ])
    HLA_HD_COLS = ['gene_id', 'Allele1', 'Allele2']
    
    alleles = []
    with open("your_file.txt") as f:
        for line in f:
            fields = line.strip().split("\t")
            gene = fields[0]
            if gene in {geneset}:
                for allele in fields[1:]:
                    if allele not in {"Not typed", "-"}:
                        alleles.append(allele)



    def create_personalized_hla_fasta(alleles_to_extract, full_hla_fasta_path, output_fasta_path):
        not_found = []
        found_alleles = set()

        try:
            
            full_hla_fasta = pysam.FastaFile(full_hla_fasta_path)
            available_references = set(full_hla_fasta.references)

            with open(output_fasta_path, 'w') as outfile:
                for allele in alleles_to_extract:
                    if allele in available_references:
                        seq = full_hla_fasta.fetch(allele)
                        outfile.write(f">{allele}\n{seq}\n")
                        found_alleles.add(allele)
                    else:
                        found = False
                        for test_allele in available_references:
                            print(test_allele,allele)
                            if allele in test_allele:
                                found =True
                                break

                        if found:
                
                            seq = full_hla_fasta.fetch(test_allele)
                            outfile.write(f">{allele}\n{seq}\n")
                            found_alleles.add(allele)
                        #If still not found maybe the allele is too long ie the allele is hla_b_45_01_01 but the ref fasta only has hla_b_45_01
                        elif len(allele)>11:
                            
                            for i in range(1,len(allele)+1):
                                try:
                                    if allele[0:-i] in available_references:
                                        found =True
                                        seq = full_hla_fasta.fetch(allele[0:-i])
                                        outfile.write(f">{allele[0:-i]}\n{seq}\n")
                                        found_alleles.add(allele)
                                        break
                                except:
                                    pass
                        else:
                            not_found.append(allele)
            full_hla_fasta.close()

        return found_alleles


    create_personalized_hla_fasta(alleles, ${hla_sequences_db}, ${prefix}_hla_reference.fasta)


    # Create versions file
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write('    python: "3.8"\\n')
    """
}