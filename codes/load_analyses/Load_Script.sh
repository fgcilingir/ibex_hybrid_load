#!/bin/bash

################################################################################
# Ibex Genetic Load Analysis Pipeline
################################################################################
#
# Description:
#   This script calculates genetic load estimates from ibex genomic data by
#   analyzing deleterious mutation burden across different Capra species.
#   It identifies:
#   - Species-specific deleterious variants
#   - Heterosis candidates in hybrids between domestic goat and Alpine ibex
#   - Loss-of-function (LOF) mutations
#
# Input Requirements:
#   - Annotated VCF file with variant effect predictions
#   - Sample ID lists for each species/population
#   - 'searchlist' file with impact categories (HIGH, MODIFIER, etc.)
#   - Expressed CDS bed files (for filtering to expressed sites)
#
# Output:
#   - VCF files filtered by impact category and ancestral state
#   - PLINK sample count files (heterozygote/homozygote counts)
#   - Site lists for derived variants for each species
#
# Dependencies:
#   - vcftools
#   - plink2
#   - tabix/bgzip
##
# Author: Christine
# Last Updated: November 2025
################################################################################

set -euo pipefail  # Exit on error, undefined variable, or pipe failure

################################################################################
# Environment Setup
################################################################################

# Load shell configuration and activate conda environment
source ~/.bashrc
source ~/.bash_profile
micromamba activate micromamba-env

# Define directory paths
work="$SCRATCH/ibex"
inds="/data/INPUT"

# Input file
annInput="ibex.ann.filt_18Feb25"

# Navigate to working directory
cd "$work/simpleLoad"

# Chromosome identifier (single file for ibex)
CHROM="ibex"

################################################################################
# Analysis 1: Domestic Goat Load Comparison
################################################################################
# Objective:
#   Identify variants that are derived (deleterious) in domestic goats but
#   ancestral in other species. This output can also be used for analyzing
#   genetic load in hybrid populations.
#
# Approach:
#   1. Calculate allele frequencies excluding domestic goats and hybrids
#   2. Identify sites fixed for reference or alternative allele
#   3. Extract variants by impact category (HIGH, MODIFIER, LOF)
#   4. Count heterozygotes and homozygotes per individual
################################################################################

# Calculate allele frequencies in all samples except domestic goats and hybrids
vcftools --gzvcf "$work/$annInput.vcf.gz" \
    --remove "$inds/dom.txt" \
    --remove "$inds/hyb.txt" \
    --freq2 \
    --stdout | gzip > WOdom.frq.gz

# Identify sites fixed for reference allele (ancestral state) in non-domestic samples
zcat WOdom.frq.gz | awk '$5 == 1 {print $1, $2}' > WOdomFixRef.pos

# Identify sites fixed for alternative allele in non-domestic samples
zcat WOdom.frq.gz | awk '$6 == 1 {print $1, $2}' > WOdomFixAlt.pos

# Create VCF subset: sites fixed for reference in non-domestic (majority of sites)
vcftools --gzvcf "$work/$annInput.vcf.gz" \
    --positions WOdomFixRef.pos \
    --recode --recode-INFO-all \
    --stdout | bgzip -c > "$CHROM.ann.WOdomFixRef.vcf.gz"
tabix "$CHROM.ann.WOdomFixRef.vcf.gz"

# Create VCF subset: sites fixed for alternative in non-domestic (minority of sites)
vcftools --gzvcf "$work/$annInput.vcf.gz" \
    --positions WOdomFixAlt.pos \
    --recode --recode-INFO-all \
    --stdout | bgzip -c > "$CHROM.ann.WOdomFixAlt.vcf.gz"
tabix "$CHROM.ann.WOdomFixAlt.vcf.gz"

# Extract VCF headers for downstream processing
zgrep "#" "$CHROM.ann.WOdomFixRef.vcf.gz" > header1a
zgrep "#" "$CHROM.ann.WOdomFixAlt.vcf.gz" > header1b

# Filter variants by impact category (from searchlist file)
# searchlist format: impact,search_term (e.g., HIGH,HIGH)
while IFS=, read -r impact search; do
    echo "Processing impact category: $impact"
    
    # Extract variants with specific impact for FixRef sites
    zgrep -v "#" "$CHROM.ann.WOdomFixRef.vcf.gz" | \
        awk -F "|" -v search="$search" '$3 == search {print $0}' > tmp
    cat header1a tmp | bgzip > "$CHROM.$impact.WOdomFixRef.vcf.gz"
    tabix "$CHROM.$impact.WOdomFixRef.vcf.gz"
    
    # Extract variants with specific impact for FixAlt sites
    zgrep -v "#" "$CHROM.ann.WOdomFixAlt.vcf.gz" | \
        awk -F "|" -v search="$search" '$3 == search {print $0}' > tmp
    cat header1b tmp | bgzip > "$CHROM.$impact.WOdomFixAlt.vcf.gz"
    tabix "$CHROM.$impact.WOdomFixAlt.vcf.gz"
done < searchlist

# Extract Loss-of-Function (LOF) variants from HIGH impact category
zgrep "LOF" "$CHROM.HIGH.WOdomFixRef.vcf.gz" | grep -v "INFO" > tmp
cat header1a tmp | bgzip > "$CHROM.LOF.WOdomFixRef.vcf.gz"
tabix "$CHROM.LOF.WOdomFixRef.vcf.gz"

zgrep "LOF" "$CHROM.HIGH.WOdomFixAlt.vcf.gz" | grep -v "INFO" > tmp
cat header1b tmp | bgzip > "$CHROM.LOF.WOdomFixAlt.vcf.gz"
tabix "$CHROM.LOF.WOdomFixAlt.vcf.gz"

# Count heterozygotes and homozygotes per individual using PLINK2
for impact in HIGH MODIFIER LOF; do
    echo "Calculating genotype counts for: $impact"
    
    # Count heterozygotes
    plink2 --vcf "$CHROM.$impact.WOdomFixAlt.vcf.gz" \
        --allow-extra-chr --sample-counts cols=het \
        --out "$CHROM.$impact.WOdomFixAlt.het"
    
    plink2 --vcf "$CHROM.$impact.WOdomFixRef.vcf.gz" \
        --allow-extra-chr --sample-counts cols=het \
        --out "$CHROM.$impact.WOdomFixRef.het"
    
    # Count homozygotes
    plink2 --vcf "$CHROM.$impact.WOdomFixAlt.vcf.gz" \
        --allow-extra-chr --sample-counts cols=homref \
        --out "$CHROM.$impact.WOdomFixAlt.homref"
    
    plink2 --vcf "$CHROM.$impact.WOdomFixRef.vcf.gz" \
        --allow-extra-chr --sample-counts cols=homalt \
        --out "$CHROM.$impact.WOdomFixRef.homalt"
done

# Filter to expressed sites only
# Combine expressed CDS bed files for filtering
cat WOdomFixAlt.positions_expressed_cds.bed \
    WOdomFixRef.positions_expressed_cds.bed | \
    awk 'OFS="\t" {print $1, $2, $3}' > WOdomFix.expressed

# Recalculate genotype counts for expressed sites only
for impact in HIGH MODIFIER LOF; do
    echo "Calculating expressed-site genotype counts for: $impact"
    
    # Note: Rare variants may be filtered out entirely
    
    # Count heterozygotes at expressed sites
    plink2 --vcf "$CHROM.$impact.WOdomFixAlt.vcf.gz" \
        --allow-extra-chr --extract bed1 WOdomFix.expressed \
        --sample-counts cols=het \
        --out "$CHROM.$impact.WOdomFixAlt.EXPR.het"
    
    plink2 --vcf "$CHROM.$impact.WOdomFixRef.vcf.gz" \
        --allow-extra-chr --extract bed1 WOdomFix.expressed \
        --sample-counts cols=het \
        --out "$CHROM.$impact.WOdomFixRef.EXPR.het"
    
    # Count homozygotes at expressed sites
    plink2 --vcf "$CHROM.$impact.WOdomFixAlt.vcf.gz" \
        --allow-extra-chr --extract bed1 WOdomFix.expressed \
        --sample-counts cols=homref \
        --out "$CHROM.$impact.WOdomFixAlt.EXPR.homref"
    
    plink2 --vcf "$CHROM.$impact.WOdomFixRef.vcf.gz" \
        --allow-extra-chr --extract bed1 WOdomFix.expressed \
        --sample-counts cols=homalt \
        --out "$CHROM.$impact.WOdomFixRef.EXPR.homalt"
done

################################################################################
# Analysis 2: Alpine Ibex (C. ibex) Load Comparison
################################################################################
# Objective:
#   Identify variants unique to Alpine ibex populations by finding sites
#   that are fixed in all other species but variable in ibex.
#
# Note: Also excludes hybrids (real hybrids only, not the "local goats")
################################################################################

# Calculate allele frequencies excluding Alpine ibex and true hybrids
vcftools --gzvcf "$work/$annInput.vcf.gz" \
    --remove "$inds/caib.txt" \
    --remove "$inds/hybWOnewGoat.txt" \
    --freq2 \
    --stdout | gzip > WOcaib.frq.gz

# Identify fixed sites in non-ibex samples
zcat WOcaib.frq.gz | awk '$5 == 1 {print $1, $2}' > WOcaibFixRef.pos
zcat WOcaib.frq.gz | awk '$6 == 1 {print $1, $2}' > WOcaibFixAlt.pos

# Create VCF subsets for sites fixed in non-ibex samples
vcftools --gzvcf "$work/$annInput.vcf.gz" \
    --positions WOcaibFixRef.pos \
    --recode --recode-INFO-all \
    --stdout | bgzip -c > "$CHROM.ann.WOcaibFixRef.vcf.gz"
tabix "$CHROM.ann.WOcaibFixRef.vcf.gz"

vcftools --gzvcf "$work/$annInput.vcf.gz" \
    --positions WOcaibFixAlt.pos \
    --recode --recode-INFO-all \
    --stdout | bgzip -c > "$CHROM.ann.WOcaibFixAlt.vcf.gz"
tabix "$CHROM.ann.WOcaibFixAlt.vcf.gz"

# Extract VCF headers
zgrep "#" "$CHROM.ann.WOcaibFixRef.vcf.gz" > header2a
zgrep "#" "$CHROM.ann.WOcaibFixAlt.vcf.gz" > header2b

# Filter by impact category
while IFS=, read -r impact search; do
    echo "Processing Alpine ibex impact category: $impact"
    
    # Extract variants with specific impact for FixRef sites
    zgrep -v "#" "$CHROM.ann.WOcaibFixRef.vcf.gz" | \
        awk -F "|" -v search="$search" '$3 == search {print $0}' > tmp2
    cat header2a tmp2 | bgzip > "$CHROM.$impact.WOcaibFixRef.vcf.gz"
    tabix "$CHROM.$impact.WOcaibFixRef.vcf.gz"
    
    # Extract variants with specific impact for FixAlt sites
    zgrep -v "#" "$CHROM.ann.WOcaibFixAlt.vcf.gz" | \
        awk -F "|" -v search="$search" '$3 == search {print $0}' > tmp2
    cat header2b tmp2 | bgzip > "$CHROM.$impact.WOcaibFixAlt.vcf.gz"
    tabix "$CHROM.$impact.WOcaibFixAlt.vcf.gz"
done < searchlist

# Extract LOF variants from HIGH impact category
zgrep "LOF" "$CHROM.HIGH.WOcaibFixRef.vcf.gz" | grep -v "INFO" > tmp
cat header2a tmp | bgzip > "$CHROM.LOF.WOcaibFixRef.vcf.gz"
tabix "$CHROM.LOF.WOcaibFixRef.vcf.gz"

zgrep "LOF" "$CHROM.HIGH.WOcaibFixAlt.vcf.gz" | grep -v "INFO" > tmp
cat header2b tmp | bgzip > "$CHROM.LOF.WOcaibFixAlt.vcf.gz"
tabix "$CHROM.LOF.WOcaibFixAlt.vcf.gz"

# Calculate genotype counts for each impact category
for impact in HIGH MODIFIER LOF; do
    echo "Calculating genotype counts for Alpine ibex: $impact"
    
    # Count heterozygotes
    plink2 --vcf "$CHROM.$impact.WOcaibFixAlt.vcf.gz" \
        --allow-extra-chr --sample-counts cols=het \
        --out "$CHROM.$impact.WOcaibFixAlt.het"
    
    plink2 --vcf "$CHROM.$impact.WOcaibFixRef.vcf.gz" \
        --allow-extra-chr --sample-counts cols=het \
        --out "$CHROM.$impact.WOcaibFixRef.het"
    
    # Count homozygotes
    plink2 --vcf "$CHROM.$impact.WOcaibFixAlt.vcf.gz" \
        --allow-extra-chr --sample-counts cols=homref \
        --out "$CHROM.$impact.WOcaibFixAlt.homref"
    
    plink2 --vcf "$CHROM.$impact.WOcaibFixRef.vcf.gz" \
        --allow-extra-chr --sample-counts cols=homalt \
        --out "$CHROM.$impact.WOcaibFixRef.homalt"
done

# Filter to expressed sites only (added November 10, 2025)
# Combine expressed CDS bed files
cat WOcaibFixAlt.positions_expressed_cds.bed \
    WOcaibFixRef.positions_expressed_cds.bed | \
    awk 'OFS="\t" {print $1, $2, $3}' > WOcaibFix.expressed

# Recalculate genotype counts for expressed sites
for impact in HIGH MODIFIER LOF; do
    echo "Calculating expressed-site genotype counts for Alpine ibex: $impact"
    
    # Count heterozygotes at expressed sites
    plink2 --vcf "$CHROM.$impact.WOcaibFixAlt.vcf.gz" \
        --allow-extra-chr --extract bed1 WOcaibFix.expressed \
        --sample-counts cols=het \
        --out "$CHROM.$impact.WOcaibFixAlt.EXPR.het"
    
    plink2 --vcf "$CHROM.$impact.WOcaibFixRef.vcf.gz" \
        --allow-extra-chr --extract bed1 WOcaibFix.expressed \
        --sample-counts cols=het \
        --out "$CHROM.$impact.WOcaibFixRef.EXPR.het"
    
    # Count homozygotes at expressed sites
    plink2 --vcf "$CHROM.$impact.WOcaibFixAlt.vcf.gz" \
        --allow-extra-chr --extract bed1 WOcaibFix.expressed \
        --sample-counts cols=homref \
        --out "$CHROM.$impact.WOcaibFixAlt.EXPR.homref"
    
    plink2 --vcf "$CHROM.$impact.WOcaibFixRef.vcf.gz" \
        --allow-extra-chr --extract bed1 WOcaibFix.expressed \
        --sample-counts cols=homalt \
        --out "$CHROM.$impact.WOcaibFixRef.EXPR.homalt"
done

################################################################################
# Analysis 3: Heterosis Candidates
################################################################################
# Objective:
#   Identify sites that are potential candidates for heterosis effects.
#   These are sites fixed for the ancestral (reference) allele in ibex but
#   derived (alternative) in all other species.
#
# Rationale:
#   Since the reference genome is an ibex, we focus on sites in WOcaibFixAlt
#   where ibex are fixed for the reference (ancestral) state while other
#   species carry the derived allele. Hybrids carrying derived alleles at
#   these sites may show heterosis.
################################################################################

# Calculate allele frequencies in Alpine ibex only
vcftools --gzvcf "$CHROM.ann.WOcaibFixAlt.vcf.gz" \
    --keep "$inds/caib.txt" \
    --freq2 \
    --stdout | gzip > WOcaib.caibfrq.gz

# Extract sites fixed for derived allele in Alpine ibex
zcat WOcaib.caibfrq.gz | awk '$5 == 1 {print $1, $2}' > CaibFixDerived.pos

# Create VCF subset for heterosis candidate sites
vcftools --gzvcf "$work/$annInput.vcf.gz" \
    --positions CaibFixDerived.pos \
    --recode --recode-INFO-all \
    --stdout | bgzip -c > "$CHROM.ann.CaibFixDerived.vcf.gz"
tabix "$CHROM.ann.CaibFixDerived.vcf.gz"

# Extract VCF header
zgrep "#" "$CHROM.ann.CaibFixDerived.vcf.gz" > header3

# Filter by impact category
while IFS=, read -r impact search; do
    echo "Processing heterosis candidate impact category: $impact"
    
    zgrep -v "#" "$CHROM.ann.CaibFixDerived.vcf.gz" | \
        awk -F "|" -v search="$search" '$3 == search {print $0}' > tmp3
    cat header3 tmp3 | bgzip > "$CHROM.$impact.CaibFixDerived.vcf.gz"
    tabix "$CHROM.$impact.CaibFixDerived.vcf.gz"
done < searchlist

# Extract LOF variants
zgrep "LOF" "$CHROM.HIGH.CaibFixDerived.vcf.gz" | grep -v "INFO" > tmp3
cat header3 tmp3 | bgzip > "$CHROM.LOF.CaibFixDerived.vcf.gz"
tabix "$CHROM.LOF.CaibFixDerived.vcf.gz"

# Calculate genotype counts (heterozygotes most important, homozygotes for validation)
for impact in HIGH MODIFIER LOF; do
    echo "Calculating genotype counts for heterosis candidates: $impact"
    
    # Count heterozygotes
    plink2 --vcf "$CHROM.$impact.CaibFixDerived.vcf.gz" \
        --allow-extra-chr --sample-counts cols=het \
        --out "$CHROM.$impact.CaibFixDerived.het"
    
    # Count homozygotes (for validation)
    plink2 --vcf "$CHROM.$impact.CaibFixDerived.vcf.gz" \
        --allow-extra-chr --sample-counts cols=homref \
        --out "$CHROM.$impact.CaibFixDerived.homref"
    
    plink2 --vcf "$CHROM.$impact.CaibFixDerived.vcf.gz" \
        --allow-extra-chr --sample-counts cols=homalt \
        --out "$CHROM.$impact.CaibFixDerived.homalt"
done

# Recalculate for expressed sites only
for impact in HIGH MODIFIER LOF; do
    echo "Calculating expressed-site counts for heterosis candidates: $impact"
    
    # Count heterozygotes at expressed sites
    plink2 --vcf "$CHROM.$impact.CaibFixDerived.vcf.gz" \
        --allow-extra-chr --extract bed1 WOcaibFix.expressed \
        --sample-counts cols=het \
        --out "$CHROM.$impact.CaibFixDerived.EXPR.het"
    
    # Count homozygotes at expressed sites (for validation)
    plink2 --vcf "$CHROM.$impact.CaibFixDerived.vcf.gz" \
        --allow-extra-chr --extract bed1 WOcaibFix.expressed \
        --sample-counts cols=homref \
        --out "$CHROM.$impact.CaibFixDerived.EXPR.homref"
    
    plink2 --vcf "$CHROM.$impact.CaibFixDerived.vcf.gz" \
        --allow-extra-chr --extract bed1 WOcaibFix.expressed \
        --sample-counts cols=homalt \
        --out "$CHROM.$impact.CaibFixDerived.EXPR.homalt"
done

################################################################################
# Summary: LOF Variant Site Lists (for Supplementary Tables)
################################################################################
# Generate site lists for Loss-of-Function variants across different categories
# for use in manuscript supplementary tables.
################################################################################

# LOF heterosis candidates (Result: 7 sites)
vcftools --gzvcf ibex.LOF.CaibFixDerived.vcf.gz \
    --kept-sites \
    --out LOF.CaibFixDerived

# LOF variants segregating among Alpine ibex (Result: 155 + 22 = 177 sites)
vcftools --gzvcf "$CHROM.LOF.WOcaibFixRef.vcf.gz" \
    --keep "$inds/caib.txt" \
    --mac 1 \
    --kept-sites \
    --out LOF.WOcaibFixRef.mac1ibex

vcftools --gzvcf "$CHROM.LOF.WOcaibFixAlt.vcf.gz" \
    --keep "$inds/caib.txt" \
    --mac 1 \
    --kept-sites \
    --out LOF.WOcaibFixAlt.mac1ibex

# Newly acquired LOF mutations in hybrids (Result: 0 sites in FixAlt)
vcftools --gzvcf ibex.LOF.WOdomFixAlt.vcf.gz

vcftools --gzvcf ibex.LOF.WOdomFixRef.vcf.gz \
    --keep "$inds/hybWOnewGoat.txt" \
    --non-ref-ac 1 \
    --kept-sites \
    --out LOF.WOdomFixRef.hybWOnewGoat.mac1

# LOF sites segregating among domestic goats (Result: 722 sites, 817 with hybrids)
vcftools --gzvcf "$CHROM.LOF.WOdomFixAlt.vcf.gz"  # No sites

vcftools --gzvcf "$CHROM.LOF.WOdomFixRef.vcf.gz" \
    --keep "$inds/dom.txt" \
    --mac 1



################################################################################
# Analysis 4: Iberian Ibex (C. pyrenaica) Load Comparison
################################################################################
# Objective:
#   Identify variants unique to C. pyrenaica (capy) populations.
################################################################################

# Calculate allele frequencies excluding C. pyrenaica
vcftools --gzvcf "$work/$annInput.vcf.gz" \
    --remove "$inds/capy.txt" \
    --freq2 \
    --stdout | gzip > WOcapy.frq.gz

# Identify fixed sites in non-capy samples
zcat WOcapy.frq.gz | awk '$5 == 1 {print $1, $2}' > WOcapyFixRef.pos
zcat WOcapy.frq.gz | awk '$6 == 1 {print $1, $2}' > WOcapyFixAlt.pos

# Create VCF subsets for fixed sites
vcftools --gzvcf "$work/$annInput.vcf.gz" \
    --positions WOcapyFixRef.pos \
    --recode --recode-INFO-all \
    --stdout | bgzip -c > "$CHROM.ann.WOcapyFixRef.vcf.gz"
tabix "$CHROM.ann.WOcapyFixRef.vcf.gz"

vcftools --gzvcf "$work/$annInput.vcf.gz" \
    --positions WOcapyFixAlt.pos \
    --recode --recode-INFO-all \
    --stdout | bgzip -c > "$CHROM.ann.WOcapyFixAlt.vcf.gz"
tabix "$CHROM.ann.WOcapyFixAlt.vcf.gz"

# Extract VCF headers
zgrep "#" "$CHROM.ann.WOcapyFixRef.vcf.gz" > header4a
zgrep "#" "$CHROM.ann.WOcapyFixAlt.vcf.gz" > header4b

# Filter by impact category
while IFS=, read -r impact search; do
    echo "Processing C. pyrenaica impact category: $impact"
    
    zgrep -v "#" "$CHROM.ann.WOcapyFixRef.vcf.gz" | \
        awk -F "|" -v search="$search" '$3 == search {print $0}' > tmp
    cat header4a tmp | bgzip > "$CHROM.$impact.WOcapyFixRef.vcf.gz"
    tabix "$CHROM.$impact.WOcapyFixRef.vcf.gz"
    
    zgrep -v "#" "$CHROM.ann.WOcapyFixAlt.vcf.gz" | \
        awk -F "|" -v search="$search" '$3 == search {print $0}' > tmp
    cat header4b tmp | bgzip > "$CHROM.$impact.WOcapyFixAlt.vcf.gz"
    tabix "$CHROM.$impact.WOcapyFixAlt.vcf.gz"
done < searchlist

# Extract LOF variants
zgrep "LOF" "$CHROM.HIGH.WOcapyFixRef.vcf.gz" | grep -v "INFO" > tmp
cat header4a tmp | bgzip > "$CHROM.LOF.WOcapyFixRef.vcf.gz"
tabix "$CHROM.LOF.WOcapyFixRef.vcf.gz"

zgrep "LOF" "$CHROM.HIGH.WOcapyFixAlt.vcf.gz" | grep -v "INFO" > tmp
cat header4b tmp | bgzip > "$CHROM.LOF.WOcapyFixAlt.vcf.gz"
tabix "$CHROM.LOF.WOcapyFixAlt.vcf.gz"

# Calculate genotype counts
for impact in HIGH MODIFIER LOF; do
    echo "Calculating genotype counts for C. pyrenaica: $impact"
    
    # Count heterozygotes
    plink2 --vcf "$CHROM.$impact.WOcapyFixAlt.vcf.gz" \
        --allow-extra-chr --sample-counts cols=het \
        --out "$CHROM.$impact.WOcapyFixAlt.het"
    
    plink2 --vcf "$CHROM.$impact.WOcapyFixRef.vcf.gz" \
        --allow-extra-chr --sample-counts cols=het \
        --out "$CHROM.$impact.WOcapyFixRef.het"
    
    # Count homozygotes
    plink2 --vcf "$CHROM.$impact.WOcapyFixAlt.vcf.gz" \
        --allow-extra-chr --sample-counts cols=homref \
        --out "$CHROM.$impact.WOcapyFixAlt.homref"
    
    plink2 --vcf "$CHROM.$impact.WOcapyFixRef.vcf.gz" \
        --allow-extra-chr --sample-counts cols=homalt \
        --out "$CHROM.$impact.WOcapyFixRef.homalt"
done

################################################################################
# Analysis 5: Siberian Ibex (C. sibirica) Load Comparison
################################################################################
# Objective:
#   Identify variants unique to C. sibirica (casi) populations.
################################################################################

# Calculate allele frequencies excluding C. sibirica
vcftools --gzvcf "$work/$annInput.vcf.gz" \
    --remove "$inds/casi.txt" \
    --freq2 \
    --stdout | gzip > WOcasi.frq.gz

# Identify fixed sites in non-casi samples
zcat WOcasi.frq.gz | awk '$5 == 1 {print $1, $2}' > WOcasiFixRef.pos
zcat WOcasi.frq.gz | awk '$6 == 1 {print $1, $2}' > WOcasiFixAlt.pos

# Create VCF subsets for fixed sites
vcftools --gzvcf "$work/$annInput.vcf.gz" \
    --positions WOcasiFixRef.pos \
    --recode --recode-INFO-all \
    --stdout | bgzip -c > "$CHROM.ann.WOcasiFixRef.vcf.gz"
tabix "$CHROM.ann.WOcasiFixRef.vcf.gz"

vcftools --gzvcf "$work/$annInput.vcf.gz" \
    --positions WOcasiFixAlt.pos \
    --recode --recode-INFO-all \
    --stdout | bgzip -c > "$CHROM.ann.WOcasiFixAlt.vcf.gz"
tabix "$CHROM.ann.WOcasiFixAlt.vcf.gz"

# Extract VCF headers
zgrep "#" "$CHROM.ann.WOcasiFixRef.vcf.gz" > header5a
zgrep "#" "$CHROM.ann.WOcasiFixAlt.vcf.gz" > header5b

# Filter by impact category
while IFS=, read -r impact search; do
    echo "Processing C. sibirica impact category: $impact"
    
    zgrep -v "#" "$CHROM.ann.WOcasiFixRef.vcf.gz" | \
        awk -F "|" -v search="$search" '$3 == search {print $0}' > tmp
    cat header5a tmp | bgzip > "$CHROM.$impact.WOcasiFixRef.vcf.gz"
    tabix "$CHROM.$impact.WOcasiFixRef.vcf.gz"
    
    zgrep -v "#" "$CHROM.ann.WOcasiFixAlt.vcf.gz" | \
        awk -F "|" -v search="$search" '$3 == search {print $0}' > tmp
    cat header5b tmp | bgzip > "$CHROM.$impact.WOcasiFixAlt.vcf.gz"
    tabix "$CHROM.$impact.WOcasiFixAlt.vcf.gz"
done < searchlist

# Extract LOF variants
zgrep "LOF" "$CHROM.HIGH.WOcasiFixRef.vcf.gz" | grep -v "INFO" > tmp
cat header5a tmp | bgzip > "$CHROM.LOF.WOcasiFixRef.vcf.gz"
tabix "$CHROM.LOF.WOcasiFixRef.vcf.gz"

zgrep "LOF" "$CHROM.HIGH.WOcasiFixAlt.vcf.gz" | grep -v "INFO" > tmp
cat header5b tmp | bgzip > "$CHROM.LOF.WOcasiFixAlt.vcf.gz"
tabix "$CHROM.LOF.WOcasiFixAlt.vcf.gz"

# Calculate genotype counts
for impact in HIGH MODIFIER LOF; do
    echo "Calculating genotype counts for C. sibirica: $impact"
    
    # Count heterozygotes
    plink2 --vcf "$CHROM.$impact.WOcasiFixAlt.vcf.gz" \
        --allow-extra-chr --sample-counts cols=het \
        --out "$CHROM.$impact.WOcasiFixAlt.het"
    
    plink2 --vcf "$CHROM.$impact.WOcasiFixRef.vcf.gz" \
        --allow-extra-chr --sample-counts cols=het \
        --out "$CHROM.$impact.WOcasiFixRef.het"
    
    # Count homozygotes
    plink2 --vcf "$CHROM.$impact.WOcasiFixAlt.vcf.gz" \
        --allow-extra-chr --sample-counts cols=homref \
        --out "$CHROM.$impact.WOcasiFixAlt.homref"
    
    plink2 --vcf "$CHROM.$impact.WOcasiFixRef.vcf.gz" \
        --allow-extra-chr --sample-counts cols=homalt \
        --out "$CHROM.$impact.WOcasiFixRef.homalt"
done