#!/bin/bash

# Function to display script usage
usage() {
  echo "Usage: $0 [OPTIONS] -n1 <name for vcf1> -n2 <name for vcf2> -f1 <vcf1> -f2 <vcf2> -o <outdir> -p <prefix>"
  echo "          -r1 [FILE]: txt file containing new sample name for vcf1 in the order of the sample name in vcf1"
  echo "          -r2 [FILE]: txt file containing new sample name for vcf2 in the order of the sample name in vcf2"
  echo "Description:"
  echo "        *  This program merges 2 vcfs with the same sample list, in favor of vcf1 when both vcfs has the same variant."
  echo "        *  This program only doing intersect, avoiding duplication, doing union while maintaining origin information."
  echo "        *  In case of sample name format not compatiable, use -r1 or r2 or both to rename samples within vcfs."
  echo "        *  Variant's origin will be indicated by adding tag in the INFO field, with -n1 or -n2, or -n1 + -n2"
  echo "        *  New tag with -n2+FILTER indicates that vcf2 contains this variant but FILTER!=PASS"
}

# Initialize variables
name1=""
name2=""
vcf1=""
vcf2=""
outdir=""
prefix=""
sample_names1=""
sample_names2=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    -n1) name1="$2"; shift; shift ;;
    -n2) name2="$2"; shift; shift ;;
    -f1) vcf1_input="$2"; shift; shift  ;;
    -f2) vcf2_input="$2"; shift; shift  ;;
    -o) outdir="$2"; shift; shift ;;
    -p) prefix="$2"; shift; shift ;;
    -r1) sample_names1="$2"; shift; shift ;;
    -r2) sample_names2="$2"; shift; shift ;;
    -h) usage; exit 0 ;;
    *) echo "Error: Invalid option $key." >&2; usage; exit 1 ;;
  esac
done

# Check for missing arguments
if [ -z "$name1" ] || [ -z "$name2" ] || [ -z "$vcf1_input" ] || [ -z "$vcf2_input" ] || [ -z "$outdir" ]; then
  echo "Error: Missing arguments." >&2
  usage
  exit 1
fi

# Check if outdir exists
if [ -d "$outdir" ]; then
  echo "Output directory '$outdir' already exists."
  read -p "Do you want to delete its contents and continue? (y/n): " confirmation
  [[ "$confirmation" != [Yy] ]] && echo "Aborted." && exit 0
  rm -rf "$outdir"/*
else
  mkdir -p "$outdir"
fi

echo

# Function to zip/unzip VCF using bgzip
prepare_vcf() {
  local vcf_file="$1"
  local outdir="$2"
  local sample_names="$3"

  if [[ "$vcf_file" == *.gz ]]; then
    local vcf_unzipped=$(basename "$vcf_file" ".gz")
    gunzip -c "$vcf_file" > "$outdir/$vcf_unzipped"
  else
    local vcf_unzipped=`basename $vcf_file`
    cp $vcf_file $outdir/$vcf_unzipped
  fi

  if [[ "$sample_names" != "" ]]; then
    bcftools reheader -s ${sample_names} -o $outdir/${vcf_unzipped}.tmp $outdir/$vcf_unzipped
    rm -rf $outdir/$vcf_unzipped
    mv $outdir/${vcf_unzipped}.tmp $outdir/$vcf_unzipped
  fi

  bgzip $outdir/$vcf_unzipped
  zipped_vcf=$outdir/${vcf_unzipped}.gz

  return_value=$zipped_vcf
  return
}

prepare_vcf "$vcf1_input" "$outdir" "$sample_names1"
vcf1=$return_value
prepare_vcf "$vcf2_input" "$outdir" "$sample_names2"
vcf2=$return_value


tabix --preset vcf $vcf1
tabix --preset vcf $vcf2

  # Get set differences of variant calls:
  # 0000: MuTect2 only
  # 0001: Strelka2 only
  # 0002: MuTect2 calls shared by Strelka2
  # 0003: Strelka2 calls shared by MuTect2
  bcftools isec \
    --output-type z \
    --prefix $outdir \
    ${vcf1} ${vcf2}


echo -e "##INFO=<ID=$name1,Number=0,Type=Flag,Description=\"Variant was called by $name1\">" > $outdir/vcf.header
echo -e "##INFO=<ID=$name2,Number=0,Type=Flag,Description=\"Variant was called by $name2\">" >> $outdir/vcf.header
echo -e "##INFO=<ID=${name2}FILTER,Number=0,Type=Flag,Description=\"Variant failed filters in ${name2}\">" >> $outdir/vcf.header

  bcftools annotate \
    --header-lines $outdir/vcf.header \
    --annotations ${outdir}/0000.vcf.gz \
    --mark-sites +${name1} \
    --output-type z \
    --output ${outdir}/0000.annot.vcf.gz \
    ${outdir}/0000.vcf.gz
  tabix --preset vcf ${outdir}/0000.annot.vcf.gz

  bcftools annotate \
    --header-lines $outdir/vcf.header \
    --annotations ${outdir}/0001.vcf.gz \
    --mark-sites +${name2} \
    --output-type z \
    --output ${outdir}/0001.annot.vcf.gz \
    ${outdir}/0001.vcf.gz
  tabix --preset vcf ${outdir}/0001.annot.vcf.gz

  bcftools annotate \
    --header-lines $outdir/vcf.header \
    --annotations ${outdir}/0002.vcf.gz \
    --mark-sites "+${name1};${name2}" \
    --output-type z \
    --output ${outdir}/0002.tmp.vcf.gz \
    ${outdir}/0002.vcf.gz
  tabix --preset vcf ${outdir}/0002.tmp.vcf.gz

  bcftools annotate \
    --annotations ${outdir}/0003.vcf.gz \
    --include 'FILTER!="PASS"' \
    --mark-sites "+${name2}FILTER" \
    -k \
    --output-type z \
    --output ${outdir}/0003.annot.vcf.gz \
    ${outdir}/0003.vcf.gz
  tabix --preset vcf ${outdir}/0003.annot.vcf.gz

  bcftools annotate \
    --annotations ${outdir}/0003.annot.vcf.gz \
    --columns +INFO,+FORMAT,${name2}FILTER \
    --output-type z \
    --output ${outdir}/0002.annot.vcf.gz \
    ${outdir}/0002.tmp.vcf.gz
  tabix --preset vcf ${outdir}/0002.annot.vcf.gz

  bcftools concat \
    --allow-overlaps \
    --rm-dups all \
    ${outdir}/0000.annot.vcf.gz \
    ${outdir}/0001.annot.vcf.gz \
    ${outdir}/0002.annot.vcf.gz | \
  bcftools sort > $outdir/$prefix.union.vcf
  bgzip $outdir/$prefix.union.vcf
  tabix --preset vcf $outdir/$prefix.union.vcf.gz

# Cleaning up:
rm -rf $outdir/000*.vcf.gz*
rm -rf $outdir/vcf.header
rm -rf ${vcf1}*
rm -rf ${vcf2}*
rm -rf $outdir/README.txt
