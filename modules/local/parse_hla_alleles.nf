process PARSE_HLA_ALLELES {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'biocontainers/python:3.11' }"

    input:
    tuple val(meta), path(hla_calls)

    output:
    tuple val(meta), path("alleles_list.txt"), emit: alleles_list
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
#!/usr/bin/env python3

geneset = set(["A", "B", "C"])

alleles = []
with open("${hla_calls}") as f:
    for line in f:
        fields = line.strip().split("\\t")
        gene = fields[0]

        if gene in geneset:
            # Process each allele separately (Allele1 and Allele2)
            for allele in fields[1:]:
                allele = allele.replace("*", "_").replace(":", "_")
                if allele not in {"Not typed", "-"} and allele not in alleles:
                    alleles.append(allele)

# Write each allele on a separate line
with open("alleles_list.txt", "w") as f:
    for allele in alleles:
        f.write(f"{allele}\\n")

print(f"Found {len(alleles)} unique alleles: {alleles}")

# Create versions file
with open("versions.yml", "w") as f:
    f.write('"${task.process}":\\n')
    f.write('    python: "3.11"\\n')
    """
}
