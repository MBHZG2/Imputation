#!/bin/bash

# Parse command-line options
while getopts p:r:o:i:t:B:b: flag; do
  case "${flag}" in
    p) plink=${OPTARG} ;;
    r) ref=${OPTARG} ;;
    o) outputdir=${OPTARG} ;;
    i) inputPath=${OPTARG} ;;
    t) target=${OPTARG} ;;
    B) beagle=${OPTARG} ;;
    b) bcftools_path=${OPTARG} ;;
  esac
done

# Function to check if a command is available
command_exists() {
  command -v "$1" >/dev/null 2>&1
}

# Check for dependencies
check_dependencies() {
  for tool in "$plink./plink" "$bcftools_path/bcftools"; do
    if ! command_exists "$tool"; then
      echo "Error: $tool not found. Please install the required tools and update the script accordingly." >&2
      exit 1
    fi
  done
}

#Check for dependencies beagle
check_dependencies() {
  if ! command_exists java; then
    echo "Error: Java not found. Please install Java to use Beagle.">&2
    exit 1
  fi

  if ! command_exists java -jar $beagle; then
    echo "Error: Beagle not found. Please install Beagle and update the script accordingly.">&2
    exit 1
  fi
}



# Function to filter panels
filter_panel() {
  local panel=$1
  local output_dir=$2
  local plink_cmd=$3

  mkdir -p "$output_dir"

  $plink_cmd./plink --bfile "$inputPath/$panel" \
    --allow-extra-chr --chr 1-18 --mind 0.1 --geno 0.1 --maf 0.01 --hwe 0.00001 --make-bed --out "$output_dir/$panel.FILTERED"
}

# Function to check duplicated markers
check_duplicated_markers() {
  local file_plink_format="$1"
  local output_dir="$2"
  local plink_cmd="$3"

  # ... (rest of the duplicated markers code) ...
}

# Redirect stdout and stderr to a log file
{
  echo "======== Start of Pipeline Log $(date) ========"
  echo "plinkPath: $plink"
  echo "reference: $ref"
  echo "outputDir: $outputdir"
  echo "inputPath: $inputPath"
  echo "target: $target"
  echo "Beagle: $beagle"
  echo "bcftoolsPath: $bcftools_path"

  # Check dependencies
  check_dependencies

  # Filter panels
  filter_panel "$ref" "$outputdir/$ref" "$plink"
  filter_panel "$target" "$outputdir/$target" "$plink"

  # Usage example for target and reference panels
  check_duplicated_markers "$target.FILTERED" "$outputdir" "$plink./plink"
  check_duplicated_markers "$ref.FILTERED" "$outputdir" "$plink./plink"

  # Work on each chromosome separately and recode as VCF. then phasing reference panels and imputation using Beagle
  for i in $(cat "$outputdir/chr"); do
    $plink./plink --bfile "$outputdir/$ref.FILTERED.NODUP" --chr "$i" --recode vcf --out "$outputdir/$ref.FILTERED.$i"
    $plink./plink --bfile "$outputdir/$target.FILTERED.NODUP" --chr "$i" --recode vcf --out "$outputdir/$target.FILTERED.$i"
    #java -jar "$beagle" nthreads=2 gt="$outputdir/$ref.FILTERED.$i.vcf" out="$outputdir/$ref.FILTERED.pahsed.$i.vcf"
    #java -jar "$beagle" nthreads=2 gt="$outputdir/$target.FILTERED.$i.vcf" ref="$outputdir/$ref.FILTERED.pahsed.$i.vcf.vcf.gz" out="$outputdir/$target.FILTERED.$i.imputed.vcf"
  done

  # Index and extract statistics for imputed VCF files
  for i in $(cat "$outputdir/chr"); do
    tabix -f -p vcf "$outputdir/$target.FILTERED.$i.imputed.vcf.vcf.gz"
  done

  for i in $(cat "$outputdir/chr"); do
    $bcftools_path/bcftools query -f \
      '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%DR2\t%AF\n' \
      "$outputdir/$target.FILTERED.$i.imputed.vcf.vcf.gz"
  done > "$outputdir/$target.stat.all.txt"

  # Generate MAF|DR2 plot and filter markers based on DR2 threshold
  /sw_results/tools/R-4.2.2/bin/Rscript \
    "$outputdir/plots.stats.bcftools.R" \
    "$outputdir/$target.stat.all.txt" 0.6 "$outputdir/$target.beagleSTAT.txt" "$outputdir/corr.plot"

  awk '{if ($8 > 0.6) print $3}' "$outputdir/$target.stat.all.txt" > "$outputdir/$target.markers.extract.txt"

  # Convert to bim, fam, bed, then merge and extract markers with DR2 > 0.6
  for i in $(cat "$outputdir/chr"); do
    $plink./plink --vcf "$outputdir/$target.FILTERED.$i.imputed.vcf.vcf.gz" \
      --make-bed --out "$outputdir/$target.$i.imputed"
  done

  for i in $(cat "$outputdir/chr"); do
    echo "$outputdir/$target.$i.imputed.bed $outputdir/$target.$i.imputed.bim $outputdir/$target.$i.imputed.fam"
  done > "$outputdir/$target.imputed.lista.merge.txt"

  $plink./plink --merge-list "$outputdir/$target.imputed.lista.merge.txt" \
    --extract "$outputdir/$target.markers.extract.txt" \
    --make-bed --out "$outputdir/$target.all.imputed"

  echo "======== End of Pipeline Log $(date) ========"
} | tee -a "$outputdir/pipeline.log" 2>&1
