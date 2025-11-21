# Per-Sample VCF Combination

## Overview

The `COMBINE_ALLELE_VCFS` module combines per-allele VCF files into a single VCF per sample, while preserving the allele-of-origin information for each variant.

## Purpose

While per-allele mutation calling provides the most accurate variant detection, you often want a consolidated view of all mutations across all HLA alleles for:
- Easier downstream analysis
- Comprehensive mutation burden assessment
- Simplified visualization
- Integration with other tools

## How It Works

### Input
- Multiple VCF files per sample (one per allele)
- Separate combinations for each variant caller:
  - Mutect2 (filtered variants)
  - Strelka SNVs
  - Strelka Indels

### Process

1. **Annotation**: Each variant is annotated with its source allele
   - Adds `HLA_ALLELE` INFO field
   - Example: `HLA_ALLELE=A_01_01_01_01`

2. **Combination**: All allele-specific VCFs are concatenated
   - Uses `bcftools concat --allow-overlaps`
   - Removes duplicate variants
   - Preserves all variant information

3. **Indexing**: Combined VCF is compressed and indexed
   - bgzip compression
   - tabix indexing

### Output

For each sample, you get **3 combined VCF files**:
1. `{sample}_mutect2.combined.vcf.gz` - All Mutect2 variants
2. `{sample}_strelka_snvs.combined.vcf.gz` - All Strelka SNV variants
3. `{sample}_strelka_indels.combined.vcf.gz` - All Strelka indel variants

Plus **statistics files**:
- `{sample}_{caller}.stats.txt` - Variant counts per allele

## VCF Format

### HLA_ALLELE INFO Field

Each variant in the combined VCF has an `HLA_ALLELE` annotation:

```vcf
##INFO=<ID=HLA_ALLELE,Number=1,Type=String,Description="HLA allele from which this variant was called">
```

Example variant:
```vcf
A_01_01_01_01    100    .    A    T    .    PASS    HLA_ALLELE=A_01_01_01_01;...
```

### Variant IDs

Variants are assigned unique IDs based on their genomic position:
```
{CHROM}_{POS}_{REF}_{ALT}
```

Example: `A_01_01_01_01_100_A_T`

## Usage Examples

### Querying Variants by Allele

Extract variants from a specific allele:
```bash
bcftools view -i 'INFO/HLA_ALLELE="A_01_01_01_01"' sample_mutect2.combined.vcf.gz
```

### Counting Variants per Allele

```bash
bcftools query -f '%INFO/HLA_ALLELE\n' sample_mutect2.combined.vcf.gz | sort | uniq -c
```

Output:
```
  15 A_01_01_01_01
  12 A_02_01_01_01
  18 B_07_02_01
  ...
```

### Getting Variant Statistics

The pipeline automatically generates statistics:
```bash
cat sample_mutect2.stats.txt
```

Output:
```
Combined VCF Statistics for sample123 (mutect2):
Total variants: 89

Variants per allele:
     15 A_01_01_01_01
     12 A_02_01_01_01
     18 B_07_02_01
     16 B_08_01_01
     14 C_01_02_01
     14 C_07_01_01
```

### Converting to Table Format

Extract key information to a table:
```bash
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/HLA_ALLELE\t%FILTER\n' \
    sample_mutect2.combined.vcf.gz > variants.tsv
```

### Filtering by Quality

Get high-quality variants only:
```bash
bcftools view -f PASS sample_mutect2.combined.vcf.gz > high_quality_variants.vcf
```

## Output Directory Structure

```
results/
├── sample123/
│   ├── combined_vcfs/
│   │   ├── sample123_mutect2.combined.vcf.gz
│   │   ├── sample123_mutect2.combined.vcf.gz.tbi
│   │   ├── sample123_strelka_snvs.combined.vcf.gz
│   │   ├── sample123_strelka_snvs.combined.vcf.gz.tbi
│   │   ├── sample123_strelka_indels.combined.vcf.gz
│   │   └── sample123_strelka_indels.combined.vcf.gz.tbi
│   ├── qc/
│   │   ├── sample123_mutect2.stats.txt
│   │   ├── sample123_strelka_snvs.stats.txt
│   │   └── sample123_strelka_indels.stats.txt
│   ├── allele_bams/
│   │   └── (per-allele BAMs)
│   └── mutect2/
│       └── (per-allele VCFs)
```

## Comparison: Per-Allele vs Combined

### Per-Allele VCFs (6 alleles)
```
sample_A_01_01_mutect2.vcf.gz    (15 variants)
sample_A_02_01_mutect2.vcf.gz    (12 variants)
sample_B_07_02_mutect2.vcf.gz    (18 variants)
sample_B_08_01_mutect2.vcf.gz    (16 variants)
sample_C_01_02_mutect2.vcf.gz    (14 variants)
sample_C_07_01_mutect2.vcf.gz    (14 variants)
```

### Combined VCF
```
sample_mutect2.combined.vcf.gz   (89 variants total)
  - Contains all variants from all alleles
  - Each variant tagged with source allele
  - Duplicates removed
```

## Duplicate Handling

If the same variant appears in multiple alleles (rare but possible):
- Only one instance is kept
- The `HLA_ALLELE` field shows the first allele where it was found
- This typically happens only with sequencing artifacts

## Integration with Downstream Tools

### VEP (Variant Effect Predictor)
```bash
vep --input_file sample_mutect2.combined.vcf.gz \
    --output_file annotated.vcf \
    --vcf \
    --custom HLA_ALLELE
```

### ANNOVAR
```bash
table_annovar.pl sample_mutect2.combined.vcf.gz humandb/ \
    -buildver hg38 \
    -protocol refGene \
    -operation g \
    -vcfinput
```

### bcftools stats
```bash
bcftools stats sample_mutect2.combined.vcf.gz > variant_stats.txt
```

## Performance Considerations

### Speed
- Combination is very fast (typically < 1 minute)
- Uses bcftools for efficient VCF manipulation
- Parallel processing not needed due to speed

### Storage
- Combined VCFs are smaller than sum of per-allele VCFs
- Compression is more efficient with more data
- Typical size: ~100-500 KB per sample

### Memory
- Low memory footprint (< 1 GB)
- Streaming processing with bcftools
- No need for high-memory nodes

## Advanced Usage

### Combining Multiple Samples

Merge combined VCFs from multiple samples:
```bash
bcftools merge \
    sample1_mutect2.combined.vcf.gz \
    sample2_mutect2.combined.vcf.gz \
    sample3_mutect2.combined.vcf.gz \
    -o cohort.vcf.gz -O z
```

### Extracting Allele-Specific Variants

Create separate VCFs for each allele from combined VCF:
```bash
for allele in A_01_01 A_02_01 B_07_02 B_08_01 C_01_02 C_07_01; do
    bcftools view -i "INFO/HLA_ALLELE=\"$allele\"" \
        sample_mutect2.combined.vcf.gz \
        -o sample_${allele}_extracted.vcf.gz -O z
done
```

### Mutation Signature Analysis

Extract mutation context for signature analysis:
```bash
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/HLA_ALLELE\n' \
    sample_mutect2.combined.vcf.gz | \
    awk '{print $1":"$2"-"$2"\t"$3">"$4"\t"$5}' \
    > mutations_for_signatures.txt
```

## Troubleshooting

### Issue: No HLA_ALLELE annotation
**Cause**: Allele name couldn't be extracted from filename
**Solution**: Check that VCF filenames follow pattern: `{sample}_{ALLELE}_{caller}.vcf.gz`

### Issue: Duplicate variants
**Cause**: Same variant called in multiple alleles
**Solution**: This is normal; bcftools removes duplicates automatically

### Issue: Missing variants
**Cause**: Variants may be in per-allele VCFs but not combined
**Solution**: Check that all per-allele VCFs are being found; check file paths

## Best Practices

1. **Keep Both**: Maintain both per-allele and combined VCFs
   - Per-allele: For detailed analysis and QC
   - Combined: For quick overview and downstream tools

2. **Use HLA_ALLELE Field**: Always filter/group by allele when analyzing
   - Mutations patterns may differ between alleles
   - Allele-specific analysis provides biological insights

3. **Verify Counts**: Check statistics files to ensure all alleles are included
   - Compare per-allele counts with combined totals
   - Look for unexpectedly low counts

4. **Document Filtering**: Record any filtering steps applied
   - Note which variants were excluded and why
   - Maintain reproducibility

## Related Documentation

- [POLYSOLVER_FILTERING.md](POLYSOLVER_FILTERING.md) - Read filtering strategies
- [Per-Allele Mutation Calling](../workflows/hlasomatic.nf) - Workflow implementation
- [bcftools documentation](http://samtools.github.io/bcftools/) - VCF manipulation

## References

- VCF Format Specification: https://samtools.github.io/hts-specs/VCFv4.3.pdf
- bcftools Manual: http://samtools.github.io/bcftools/bcftools.html
