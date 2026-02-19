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

geneset = set(["A", "B", "C"])
HLA_HD_COLS = ['gene_id', 'Allele1', 'Allele2']

alleles = []
with open("${hla_calls}") as f:
    for line in f:
        fields = line.strip().split("\\t")
        gene = fields[0]

        if gene in geneset:
            for allele in fields[1:]:
                allele = allele.replace("*", "_").replace(":", "_")
                if allele not in {"Not typed", "-"} and allele not in alleles:
                    alleles.append(allele)


print(alleles)
def create_personalized_hla_fasta(alleles_to_extract, full_hla_fasta_path, output_fasta_path):
    not_found = []
    found_alleles = set()

    try:

        full_hla_fasta = pysam.FastaFile(full_hla_fasta_path)
        available_references = set(full_hla_fasta.references)

        with open(output_fasta_path, 'w') as outfile:
            for allele in alleles_to_extract:
                allele = allele[4:]
                print(allele)
                if allele in available_references:
                    found = True
                    seq = full_hla_fasta.fetch(allele)
                    outfile.write(f">{allele}\\n{seq}\\n")
                    found_alleles.add(allele)
                    print(f"1 found {allele}")
                else:
                    print("???")
                    found = False
                    for test_allele in available_references:
                        if allele in test_allele:
                            found = True
                            seq = full_hla_fasta.fetch(test_allele)
                            outfile.write(f">{test_allele}\\n{seq}\\n")
                            found_alleles.add(allele)
                            print(f"2 found {allele}")
                            break

                    if len(allele) > 11 and not found:

                        for i in range(1, len(allele) + 1):
                            print(allele[0:-i])
                            try:
                                if allele[0:-i] in available_references:
                                    found = True
                                    seq = full_hla_fasta.fetch(allele[0:-i])
                                    outfile.write(f">{allele[0:-i]}\\n{seq}\\n")
                                    found_alleles.add(allele)
                                    print(f"3 found {allele}")
                                    break
                            except:
                                pass

                    if not found:
                        print(f"not found {allele}")
                        not_found.append(allele)
        full_hla_fasta.close()
        print(not_found)
    except Exception as e:
        print(e)
        pass
    return found_alleles


def append_imgt_decoys(output_fasta_path, imgt_fasta_path, typed_genes):
    \"\"\"
    Append one representative decoy allele per non-typed HLA gene from the
    IMGT/HLA reference. This covers ALL HLA genes (Class I non-classical,
    Class II, MIC, pseudogenes, etc.) as read sinks to prevent cross-gene
    contamination of A/B/C alignments.

    For each gene not in the patient's typed set, picks the first *01:01
    allele (or first available) as the representative decoy.
    \"\"\"
    imgt = pysam.FastaFile(imgt_fasta_path)
    all_refs = list(imgt.references)

    # Group alleles by gene (first field before '_')
    from collections import defaultdict
    gene_alleles = defaultdict(list)
    for ref in all_refs:
        gene = ref.split("_")[0]
        gene_alleles[gene].append(ref)

    typed_upper = set(g.upper() for g in typed_genes)
    decoys_added = 0

    with open(output_fasta_path, 'a') as outfile:
        for gene in sorted(gene_alleles.keys()):
            if gene.upper() in typed_upper:
                continue

            candidates = gene_alleles[gene]
            # Prefer *01:01 (i.e. gene_01_01) allele as representative
            chosen = None
            for allele in candidates:
                # Match gene_01_01 prefix (any trailing fields ok)
                if allele.startswith(f"{gene}_01_01"):
                    chosen = allele
                    break
            if chosen is None:
                chosen = candidates[0]

            seq = imgt.fetch(chosen)
            decoy_name = f"{chosen}_decoy"
            outfile.write(f">{decoy_name}\\n{seq}\\n")
            decoys_added += 1
            print(f"Added decoy: {decoy_name} ({len(seq)} bp)")

    imgt.close()
    print(f"Added {decoys_added} IMGT-based decoy sequences to reference")


# Step 1: Create personalized reference with patient's typed A/B/C alleles
create_personalized_hla_fasta(alleles, "${reference_fasta}", "${prefix}_hla_reference.fasta")

# Step 2: Append IMGT-based decoy sequences for all non-typed HLA genes
typed_genes = set(["A", "B", "C"])
append_imgt_decoys("${prefix}_hla_reference.fasta", "${reference_fasta}", typed_genes)


# Create versions file
with open("versions.yml", "w") as f:
    f.write('"${task.process}":\\n')
    f.write('    python: "3.8"\\n')
    """
}