process FILTER_ALLELE_BAM {
    tag "${meta.id}_${meta.allele_safe}"
    label 'process_medium'

    conda "bioconda::samtools=1.19.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0' :
        'biocontainers/samtools:1.19.2--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${prefix}.filtered.bam"), path("${prefix}.filtered.bam.bai"), emit: bam
    path "${prefix}.filter_stats.txt", emit: stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def max_events = task.ext.max_events ?: 10  // Default: max 10 mismatches+indels
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    # POLYSOLVER-style event count filtering
    # Counts mismatches (NM tag) + insertions + deletions from CIGAR
    # Filters reads exceeding max_events threshold

    samtools view -h ${bam} | \\
        awk -v max_events=${max_events} '
        BEGIN {
            OFS="\\t"
            total_reads=0
            filtered_reads=0
            removed_reads=0
        }
        /^@/ {
            print
            next
        }
        {
            total_reads++

            # Extract NM tag (edit distance = mismatches)
            nm = 0
            for (i=12; i<=NF; i++) {
                if (\$i ~ /^NM:i:/) {
                    split(\$i, a, ":")
                    nm = a[3]
                    break
                }
            }

            # Parse CIGAR string to count insertions and deletions
            cigar = \$6
            indels = 0

            # Count insertions (I) and deletions (D)
            while (match(cigar, /([0-9]+)([ID])/)) {
                len = substr(cigar, RSTART, RLENGTH-1)
                op = substr(cigar, RSTART+RLENGTH-1, 1)
                indels += len
                cigar = substr(cigar, RSTART+RLENGTH)
            }

            # Total events = mismatches + indels
            events = nm + indels

            # Filter based on event count
            if (events <= max_events) {
                # Set mapping quality to 70 (POLYSOLVER convention)
                \$5 = 70
                print
                filtered_reads++
            } else {
                removed_reads++
            }
        }
        END {
            print "Total reads:", total_reads > "/dev/stderr"
            print "Reads passing filter:", filtered_reads > "/dev/stderr"
            print "Reads removed:", removed_reads > "/dev/stderr"
            print "Filter rate:", (removed_reads/total_reads)*100 "%" > "/dev/stderr"
        }
        ' | \\
        samtools view -b -o ${prefix}.filtered.bam -

    # Generate filter statistics
    cat > ${prefix}.filter_stats.txt <<EOF
Sample: ${meta.id}
Allele: ${meta.allele}
Max events threshold: ${max_events}
Statistics written during filtering (see log)
EOF

    # Index the filtered BAM
    samtools index ${prefix}.filtered.bam

    cat <<END_VERSIONS > versions.yml
"${task.process}":
    samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    awk: \$(awk --version 2>&1 | head -n1 | sed 's/GNU Awk //' | sed 's/,.*//')
END_VERSIONS
    """
}
