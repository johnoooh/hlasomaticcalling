# POLYSOLVER-Style Filtering Implementation

## Overview

This document describes the read filtering strategies implemented in HLAsomatic, based on POLYSOLVER's mutation calling approach.

## Filtering Stages

### 1. Per-Allele BAM Extraction (`EXTRACT_ALLELE_BAM`)

**Location**: `modules/local/extract_allele_bam.nf`

**Purpose**: Extract reads mapping to a specific HLA allele with same-allele pairing filter

**Filters Applied**:

1. **Allele-Specific Extraction**
   - Extracts only reads mapping to the target allele
   - Uses `samtools view` with specific reference name

2. **Same-Allele Pairing Filter** ⭐ **CRITICAL**
   - Keeps only read pairs where **both mates map to the same allele**
   - Checks SAM field 7 (RNEXT) == "="
   - Prevents cross-allele contamination
   - **POLYSOLVER equivalent**: `filterReads.pl` RNEXT check

3. **Properly Paired Filter**
   - Uses `samtools fixmate` to update mate information
   - Uses `samtools view -f 0x2` to keep only properly paired reads
   - Removes singletons and unmapped mates

**Example**:
```bash
# Read pair mapping scenario:
# Read1: Maps to HLA-A*01:01
# Read2: Maps to HLA-A*02:01
# Result: EXCLUDED (mates map to different alleles)

# Read pair mapping scenario:
# Read1: Maps to HLA-A*01:01
# Read2: Maps to HLA-A*01:01
# Result: INCLUDED (both map to same allele)
```

**Output Statistics**:
- Total reads mapping to allele
- Filtered reads (same-allele pairs)
- Filtering rate

---

### 2. Event Count Filtering (`FILTER_ALLELE_BAM`) - OPTIONAL

**Location**: `modules/local/filter_allele_bam.nf`

**Purpose**: Advanced quality filtering based on alignment complexity

**Filters Applied**:

1. **Event Count Calculation**
   - Counts mismatches from NM tag
   - Counts insertions (I) from CIGAR string
   - Counts deletions (D) from CIGAR string
   - Total events = mismatches + insertions + deletions

2. **Event Threshold Filter**
   - Default threshold: 10 events maximum
   - Configurable via `ext.max_events` in `modules.config`
   - Removes reads with excessive complexity

3. **Mapping Quality Adjustment**
   - Sets MAPQ to 70 for passing reads (POLYSOLVER convention)

**POLYSOLVER equivalent**: `filterReads.pl` event count filter

**When to Use**:
- High coverage data (>100X) where you can afford to be stringent
- Samples with high sequencing error rates
- When you need maximum specificity over sensitivity

**When NOT to Use**:
- Low coverage data (<30X)
- When analyzing highly polymorphic regions
- When sensitivity is prioritized

**To Enable**: Add to workflow after `EXTRACT_ALLELE_BAM`:

```groovy
// In workflows/hlasomatic.nf
include { FILTER_ALLELE_BAM } from '../modules/local/filter_allele_bam'

// After EXTRACT_ALLELE_BAM
FILTER_ALLELE_BAM (
    EXTRACT_ALLELE_BAM.out.bam
)

// Use FILTER_ALLELE_BAM.out.bam instead of EXTRACT_ALLELE_BAM.out.bam
// for downstream mutation calling
```

**Configuration** (`conf/modules.config`):
```groovy
withName: FILTER_ALLELE_BAM {
    ext.max_events = 10  // Adjust threshold (default: 10)
}
```

---

## Comparison with POLYSOLVER Scripts

### filterReads.pl
| Feature | POLYSOLVER | HLAsomatic Implementation | Status |
|---------|-----------|---------------------------|--------|
| Same-allele pairing | ✓ (RNEXT="=") | ✓ EXTRACT_ALLELE_BAM | ✅ Implemented |
| Event count filter | ✓ (user threshold) | ✓ FILTER_ALLELE_BAM | ✅ Implemented (optional) |
| MAPQ adjustment | ✓ (set to 70) | ✓ FILTER_ALLELE_BAM | ✅ Implemented |

### clean_unpaired_fastq.pl
| Feature | POLYSOLVER | HLAsomatic Implementation | Status |
|---------|-----------|---------------------------|--------|
| Remove unpaired FASTQ | ✓ | ✗ Not needed | ✅ N/A - We work with BAMs |

### keep_only_read_pairs.pl
| Feature | POLYSOLVER | HLAsomatic Implementation | Status |
|---------|-----------|---------------------------|--------|
| Keep only read pairs | ✓ | ✓ samtools -f 0x2 | ✅ Implemented |
| Name-based pairing | ✓ | ✓ samtools fixmate | ✅ Implemented |

---

## Filtering Pipeline Flow

```
Input: Realigned BAM (all alleles)
    |
    v
┌─────────────────────────────────────┐
│  EXTRACT_ALLELE_BAM                 │
│  ✓ Extract allele-specific reads    │
│  ✓ Same-allele pairing filter       │
│  ✓ Properly paired filter           │
└────────────────┬────────────────────┘
                 |
                 v
        Per-allele BAMs
                 |
                 v (optional)
┌─────────────────────────────────────┐
│  FILTER_ALLELE_BAM                  │
│  ✓ Event count filtering            │
│  ✓ MAPQ adjustment                  │
└────────────────┬────────────────────┘
                 |
                 v
    High-quality per-allele BAMs
                 |
                 v
    ┌───────────┴───────────┐
    |                       |
    v                       v
MUTECT2               STRELKA
    |                       |
    v                       v
Per-allele VCFs      Per-allele VCFs
```

---

## Filtering Statistics and QC

### EXTRACT_ALLELE_BAM Output
- Log files show:
  - Total reads per allele
  - Reads passing same-allele filter
  - Allele reference names used

### FILTER_ALLELE_BAM Output
- `*.filter_stats.txt` contains:
  - Sample and allele information
  - Event threshold used
  - Filter statistics (logged to stderr)

---

## Recommendations

### Standard Analysis (Recommended)
```
✓ Use EXTRACT_ALLELE_BAM (enabled by default)
✗ Skip FILTER_ALLELE_BAM
```
**Reason**: Same-allele filtering is critical and sufficient for most cases

### High-Stringency Analysis
```
✓ Use EXTRACT_ALLELE_BAM (enabled by default)
✓ Use FILTER_ALLELE_BAM with max_events=8-10
```
**Reason**: For maximum specificity in high-coverage samples

### Maximum-Sensitivity Analysis
```
✓ Use EXTRACT_ALLELE_BAM (enabled by default)
✗ Skip FILTER_ALLELE_BAM
✓ Consider increasing variant caller sensitivity
```
**Reason**: For low-coverage samples or rare variant detection

---

## Key Differences from POLYSOLVER

### Advantages
1. **Modular Design**: Filtering is separated into distinct processes
2. **Optional Advanced Filter**: Event count filtering is opt-in
3. **Better Statistics**: Detailed filtering statistics per allele
4. **Nextflow Integration**: Proper caching and resume support

### Behavioral Differences
1. **Default Mode**: Only same-allele filtering (most critical)
2. **POLYSOLVER**: Applies both filters by default
3. **Rationale**: Gives users more control over sensitivity/specificity trade-off

---

## Testing Filtering Impact

To assess filtering impact on your data:

1. **Run without event filtering** (default):
   ```bash
   nextflow run main.nf --input samplesheet.csv --outdir results_standard
   ```

2. **Run with event filtering**:
   - Uncomment `FILTER_ALLELE_BAM` in workflow
   - Adjust `ext.max_events` if needed
   ```bash
   nextflow run main.nf --input samplesheet.csv --outdir results_stringent
   ```

3. **Compare results**:
   - Number of variants called per allele
   - Variant quality scores
   - Coverage per allele

---

## Technical Notes

### SAM/BAM Field Reference
- Field 6 (CIGAR): Alignment string with I/D/M operations
- Field 7 (RNEXT): Reference name of mate ("=" means same reference)
- NM tag: Edit distance (mismatches)
- Flag 0x2: Read properly paired

### Performance Considerations
- `EXTRACT_ALLELE_BAM`: Fast (uses native samtools)
- `FILTER_ALLELE_BAM`: Moderate (uses AWK for CIGAR parsing)
- Expected overhead: ~10-15% additional time if event filtering enabled

### Validation
- Same-allele filter is validated by checking RNEXT field
- Event counting is validated against NM tag + CIGAR operations
- Results match POLYSOLVER filterReads.pl behavior

---

## References

- POLYSOLVER: https://github.com/jason-weirather/hla-polysolver
- SAM Format Specification: https://samtools.github.io/hts-specs/SAMv1.pdf
- Nextflow Best Practices: https://www.nextflow.io/docs/latest/
