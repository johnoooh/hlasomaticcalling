process HLA_FROM_SAMPLESHEET {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'biocontainers/python:3.11' }"

    input:
    tuple val(meta), val(A1), val(A2), val(B1), val(B2), val(C1), val(C2)

    output:
    tuple val(meta), path("${prefix}_final.result.txt"), emit: hla_calls
    path "versions.yml",                                 emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
#!/usr/bin/env python3

def normalize_allele(allele, gene):
    a = allele.strip()
    if not a or a in {'-', 'Not typed', 'NA', ''}:
        return None
    if a.startswith('HLA-'):
        return a
    if a.startswith(gene + '*'):
        return 'HLA-' + a
    # bare digits e.g. "02:01:01" -> "HLA-A*02:01:01"
    return f'HLA-{gene}*{a}'

genes = [
    ('A', '${A1}', '${A2}'),
    ('B', '${B1}', '${B2}'),
    ('C', '${C1}', '${C2}'),
]

rows = []
for gene, raw1, raw2 in genes:
    allele1 = normalize_allele(raw1, gene) or 'Not typed'
    allele2 = normalize_allele(raw2, gene) or 'Not typed'
    rows.append(f'{gene}\\t{allele1}\\t{allele2}')

with open('${prefix}_final.result.txt', 'w') as f:
    f.write('\\n'.join(rows) + '\\n')

import sys
python_version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
with open('versions.yml', 'w') as f:
    f.write('"${task.process}":\\n')
    f.write(f'    python: "{python_version}"\\n')
    """
}
