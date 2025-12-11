# Ibex Genetic Load Analysis Pipeline

## Overview

This pipeline calculates genetic load estimates from ibex genomic data by analyzing deleterious mutation burden across different *Capra* species. The analysis identifies:

- Species-specific deleterious variants
- Heterosis candidates in hybrid populations  
- Loss-of-function (LOF) mutations across populations

## Requirements

### Software Dependencies

- **vcftools** - VCF file manipulation and filtering
- **plink2** - Genotype analysis and sample counting
- **tabix/bgzip** - VCF compression and indexing
- Standard Unix tools: awk, grep, zcat

### Input Files

1. **Annotated VCF file** (`ibex.ann.filt_18Feb25.vcf.gz`)
   - Variant effect predictions from annotation software (e.g., SnpEff)
   - Should include impact categories (HIGH, MODIFIER, etc.)

2. **Sample ID lists** (in `$inds` directory):
   - `dom.txt` - Domestic goats
   - `caib.txt` - Alpine ibex (*C. ibex*)
   - `hyb.txt` - All hybrids
   - `hybWOnewGoat.txt` - True hybrids only
   - `capy.txt` - Iberian ibex (*C. pyrenaica*)
   - `casi.txt` - Siberian ibex (*C. sibirica*)

3. **Impact categories file** (`searchlist`)
   - CSV format: `impact,search_term`
   - Example: `HIGH,HIGH`

4. **Expressed CDS bed files** (optional, for filtering):
   - `WOdomFix*.positions_expressed_cds.bed`
   - `WOcaibFix*.positions_expressed_cds.bed`

## Pipeline Workflow

### 1. Domestic Goat Analysis

Identifies variants that are derived (deleterious) in domestic goats but ancestral in other species.

**Outputs:**
- `ibex.*.WOdomFixRef.vcf.gz` - Sites fixed for reference in non-domestic samples
- `ibex.*.WOdomFixAlt.vcf.gz` - Sites fixed for alternative in non-domestic samples
- `*.het` and `*.homref`/`*.homalt` - Genotype counts per individual

### 2. Alpine Ibex Analysis

Identifies variants unique to Alpine ibex populations.

**Outputs:**
- `ibex.*.WOcaibFixRef.vcf.gz`
- `ibex.*.WOcaibFixAlt.vcf.gz`
- `*.het` and `*.homref`/`*.homalt` - Genotype counts per individual

### 3. Heterosis Candidate Analysis

Identifies sites potentially contributing to heterosis in hybrid populations.

**Outputs:**
- `ibex.*.CaibFixDerived.vcf.gz`
- LOF heterosis candidate list (7 sites identified)

### 4. Iberian Ibex Analysis

Identifies variants unique to *C. pyrenaica* populations.

**Outputs:**
- `ibex.*.WOcapyFixRef.vcf.gz`
- `ibex.*.WOcapyFixAlt.vcf.gz`
- `*.het` and `*.homref`/`*.homalt` - Genotype counts per individual

### 5. Siberian Ibex Analysis

Identifies variants unique to *C. sibirica* populations.

**Outputs:**
- `ibex.*.WOcasiFixRef.vcf.gz`
- `ibex.*.WOcasiFixAlt.vcf.gz`
- `*.het` and `*.homref`/`*.homalt` - Genotype counts per individual

## Output Files

### VCF Files

Impact categories:
- `*.HIGH.*` - High impact variants
- `*.MODIFIER.*` - Modifier variants
- `*.LOF.*` - Loss-of-function variants

### PLINK Sample Count Files

- `*.het.scount` - Heterozygote counts per individual
- `*.homref.scount` - Homozygous reference counts
- `*.homalt.scount` - Homozygous alternative counts
- `*.EXPR.*` - Filtered to expressed sites only

### Site Lists

- `*.kept.sites` - Lists of variant positions passing filters

## Notes

- The reference genome is an ibex, which affects interpretation of ancestral vs. derived states
- Some rare variant categories may have no variants after expression filtering

## Citation

If using this pipeline, please cite: Çilingir et al. 2025

## Contact

Christine.Grossen@wsl.ch
Last Updated: November 2025
