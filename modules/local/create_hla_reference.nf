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
    path genome_fasta
    path genome_fasta_fai

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


def append_decoy_sequences(output_fasta_path, genome_fasta_path):
    \"\"\"
    Append non-classical HLA gene and pseudogene sequences from the genome
    reference as decoy contigs. These act as read sinks to prevent paralog
    reads (HLA-E, F, G, H, J, K, L) from contaminating classical allele
    alignments and producing false positive somatic mutation calls.
    \"\"\"
    # Decoy regions for GRCh37/hg19 (no chr prefix)
    decoy_regions_grch37 = {
        "HLA_E_decoy": ("6", 30457244, 30461982),
        "HLA_F_decoy": ("6", 29690552, 29706305),
        "HLA_G_decoy": ("6", 29794744, 29798902),
        "HLA_H_decoy": ("6", 29855350, 29858259),
        "HLA_J_decoy": ("6", 29974360, 29977733),
        "HLA_K_decoy": ("6", 29894236, 29897009),
        "HLA_L_decoy": ("6", 30227339, 30234728),
    }

    # Decoy regions for GRCh38/hg38 (chr prefix)
    decoy_regions_hg38 = {
        "HLA_E_decoy": ("chr6", 30489503, 30494205),
        "HLA_F_decoy": ("chr6", 29722738, 29738528),
        "HLA_G_decoy": ("chr6", 29826967, 29831125),
        "HLA_H_decoy": ("chr6", 29887752, 29890482),
        "HLA_J_decoy": ("chr6", 30006606, 30009539),
        "HLA_K_decoy": ("chr6", 29926459, 29929232),
        "HLA_L_decoy": ("chr6", 30259562, 30266951),
    }

    genome = pysam.FastaFile(genome_fasta_path)
    chroms = set(genome.references)

    # Auto-detect genome build from chromosome naming convention
    if "chr6" in chroms:
        build = "hg38"
        decoy_regions = decoy_regions_hg38
    elif "6" in chroms:
        build = "grch37"
        decoy_regions = decoy_regions_grch37
    else:
        print("WARNING: Cannot detect genome build (neither 'chr6' nor '6' found). Skipping decoys.")
        genome.close()
        return

    print(f"Detected genome build: {build}")
    decoys_added = 0

    with open(output_fasta_path, 'a') as outfile:
        for name, (chrom, start, end) in decoy_regions.items():
            try:
                seq = genome.fetch(chrom, start, end)
                if len(seq) > 0:
                    outfile.write(f">{name}\\n{seq}\\n")
                    decoys_added += 1
                    print(f"Added decoy: {name} ({chrom}:{start}-{end}, {len(seq)} bp)")
                else:
                    print(f"WARNING: Empty sequence for {name} ({chrom}:{start}-{end})")
            except Exception as e:
                print(f"WARNING: Could not extract {name} ({chrom}:{start}-{end}): {e}")

    genome.close()
    print(f"Added {decoys_added} decoy sequences to reference")


def append_hla_y_decoy(output_fasta_path):
    \"\"\"
    HLA-Y is absent from GRCh37/GRCh38 genome assemblies (71-87% of haplotypes
    carry the deletion). Add HLA-Y CDS from IMGT/HLA as an explicit decoy to
    prevent HLA-Y reads from contaminating HLA-A alignments, especially at
    CDS position 144 where Y differs from A.
    Source: IMGT/HLA Y*01:01:01:01 (HLA:HLA13320), 1098 bp
    \"\"\"
    hla_y_cds = (
        "ATGGCGGTCGTGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGCCCTGACC"
        "CAGACCTGGGCGGGCTCCCACTCCATGAGGTATTTCTCCACATCCGTGTCCCGGCCCGGC"
        "AGTGGAGAGCCCCGCTTCATCGCAGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTC"
        "GACAGCGACGCCGCGAGCCAGAGGATGGAGCCGCGGGCGCCGTGGATGGAGCAGGAGGAG"
        "CCGGAGTATTGGGACCGGCAGACACAGATCTCCAAGACCAACGCACAGATTGACCTAGAG"
        "AGCCTGCGGATCGCGCTCCGCTACTACAACCAGAGCGAGGCCGGTTCTCACACCATCCAG"
        "AGGATGTCTGGCTGCGACGTGGGGTCGGACGGGCGCTTCCTCCGCGGGTACCGGCAGGAC"
        "GCCTACGACGGCAAGGATTACATCGCCCTGAACGAGGACCTGCGCTCTTGGACCGCGGCG"
        "GACATGGCGGCTCAGATCACCCAGCGCAAGTGGGAGGCGGCCCGTCAGGCGGAGCAGTTG"
        "AGAGCCTACCTGGAGGGCGAGTGCATGGAGTGGCTCCGCAGATACCTGGAGAACGGGAAG"
        "GAGACGCTGCAGCGCACGGACGCCCCCAAGACGCATATGACTCACCACGCTGTCTCTGAC"
        "AATGAGGCCACCCTGAGGTGCTGAGCCCTGAGCTTCTACCCTGCGGAGATCACACTGACC"
        "TGGCAGCGGGATGGGGAGGACCAGACCCAGGACACGGAGCTCGTGGAGACCAGGCCTGCA"
        "GGGGATGGAATCTTCCAGAAGTGGGCGGCTGTGGTGGTGCCTTCTGGAGAGGAGCAGAGA"
        "TACACCTGCCATGTGCAGCATGAGGGTCTGCCCAAGCCCCTCACCCTGAGATGGGAGCCG"
        "TCTTCCCATCCCACCATCCCCATCGTGGGCATCCTTGCTGGCCTGGTTCTCTTTGGAGCT"
        "GTGATCGCTGGAGCTGTGGTCGCTGCTGTGATGTGGAGGAGGAAGAGCTCAGATAGAAAA"
        "GGAGGGAGCTACTCTCAGGCTGCAAGCAGTGACATTGCCCAGGGCTCTGATGTGTCTCTC"
        "ACAGCTTGTAAAGTGTGA"
    )
    with open(output_fasta_path, 'a') as f:
        f.write(f">HLA_Y_decoy\\n{hla_y_cds}\\n")
    print(f"Added HLA-Y decoy (Y*01:01, {len(hla_y_cds)} bp) - not in genome assembly")


# Step 1: Create personalized reference with patient's typed A/B/C alleles
create_personalized_hla_fasta(alleles, "${reference_fasta}", "${prefix}_hla_reference.fasta")

# Step 2: Append non-classical HLA decoy sequences from genome reference
append_decoy_sequences("${prefix}_hla_reference.fasta", "${genome_fasta}")

# Step 3: Append HLA-Y decoy (absent from genome assembly, sourced from IMGT/HLA)
append_hla_y_decoy("${prefix}_hla_reference.fasta")


# Create versions file
with open("versions.yml", "w") as f:
    f.write('"${task.process}":\\n')
    f.write('    python: "3.8"\\n')
    """
}