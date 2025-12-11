# ============================================================================
# Analysis of Loss-of-Function (LoF) variants in Capra species
# ============================================================================
#
# Purpose: Analyze expressed LoF variants in ibex and domestic goat genomes
#          to quantify genetic load differences between species
#
# Date: November 10, 2025
# Author: Christine Grossen
#
# Input files: sampleList_v2_mod.csv and .scount files (variant counts, as produced by Load_Script.sh)
# Output: Publication-quality figures (PDF) and statistical tests
#
# ============================================================================

# Load required libraries
library(ggplot2)
library(stats)

# ============================================================================
# Configuration: Set data directory and load sample metadata
# ============================================================================
# Users should set this to the directory containing the input data files
DATA_DIR <- "./data"  # Update this path as needed

# Load sample metadata
# Expected columns: Accession, Phenotype, Category, IndID, PropGoat
info <- read.csv(file.path(DATA_DIR, "sampleList_v2_mod.csv"))

# ============================================================================
# Analysis 1: Load among domestic goat ("goat-load")
# ============================================================================
dat <- read.table(file.path(DATA_DIR, "ibex.HIGH.WOdomFixRef.EXPR.het.scount"))
dat$HighHet <- dat$V2
dat$HighHomDerived <- read.table(file.path(DATA_DIR, "ibex.HIGH.WOdomFixRef.EXPR.homalt.scount"))$V2

dat$LofHet <- read.table(file.path(DATA_DIR, "ibex.LOF.WOdomFixRef.EXPR.het.scount"))$V2
dat$LofHomDerived <- read.table(file.path(DATA_DIR, "ibex.LOF.WOdomFixRef.EXPR.homalt.scount"))$V2

dat$ModifHet <- read.table(file.path(DATA_DIR, "ibex.MODIFIER.WOdomFixRef.EXPR.het.scount"))$V2
dat$ModifHomDerived <- read.table(file.path(DATA_DIR, "ibex.MODIFIER.WOdomFixRef.EXPR.homalt.scount"))$V2

dat$HighTot <- dat$HighHet + 2 * dat$HighHomDerived
dat$LofTot <- dat$LofHet + 2 * dat$LofHomDerived
dat$ModifTot <- dat$ModifHet + 2 * dat$ModifHomDerived

dati <- merge(dat, info, by.x = "V1", by.y = "Accession")
dati$Category <- factor(dati$Category, 
                        levels = c("Cahi", "CahiNew", "hybrids", "Caib", "Capy", "other"))

datigoat <- dati
datigoat$filt <- "WOdomFixRef"

# ============================================================================
# Analysis 2: Load among Alpine ibex ("ibex-load")
# ============================================================================
dat <- read.table(file.path(DATA_DIR, "ibex.HIGH.WOcaibFixRef.EXPR.het.scount"))
dat$HighHet <- dat$V2 + read.table(file.path(DATA_DIR, "ibex.HIGH.WOcaibFixAlt.EXPR.het.scount"))$V2
dat$HighHomDerived <- read.table(file.path(DATA_DIR, "ibex.HIGH.WOcaibFixRef.EXPR.homalt.scount"))$V2 + 
                      read.table(file.path(DATA_DIR, "ibex.HIGH.WOcaibFixAlt.EXPR.homref.scount"))$V2

dat$LofHet <- read.table(file.path(DATA_DIR, "ibex.LOF.WOcaibFixRef.EXPR.het.scount"))$V2 + 
              read.table(file.path(DATA_DIR, "ibex.LOF.WOcaibFixAlt.EXPR.het.scount"))$V2
dat$LofHomDerived <- read.table(file.path(DATA_DIR, "ibex.LOF.WOcaibFixRef.EXPR.homalt.scount"))$V2 + 
                     read.table(file.path(DATA_DIR, "ibex.LOF.WOcaibFixAlt.EXPR.homref.scount"))$V2

dat$ModifHet <- read.table(file.path(DATA_DIR, "ibex.MODIFIER.WOcaibFixRef.EXPR.het.scount"))$V2 + 
                read.table(file.path(DATA_DIR, "ibex.MODIFIER.WOcaibFixAlt.EXPR.het.scount"))$V2
dat$ModifHomDerived <- read.table(file.path(DATA_DIR, "ibex.MODIFIER.WOcaibFixRef.EXPR.homalt.scount"))$V2 + 
                       read.table(file.path(DATA_DIR, "ibex.MODIFIER.WOcaibFixAlt.EXPR.homref.scount"))$V2

dat$HighTot <- dat$HighHet + 2 * dat$HighHomDerived
dat$LofTot <- dat$LofHet + 2 * dat$LofHomDerived
dat$ModifTot <- dat$ModifHet + 2 * dat$ModifHomDerived

dati <- merge(dat, info, by.x = "V1", by.y = "Accession")
dati$Category <- factor(dati$Category, 
                        levels = c("Cahi", "CahiNew", "hybrids", "Caib", "Capy", "other"))

datiibex <- dati
datiibex$filt <- "WOcaibFixRef"

# ============================================================================
# Combined Analysis: Standardized LoF load across both reference states
# ============================================================================
datboth <- rbind(datigoat, datiibex)
datboth$LofStd <- datboth$LofTot / datboth$ModifTot 

# ============================================================================
# Statistical Test: Compare standardized LoF between domestic goat and ibex
# ============================================================================

lof_test_result <- t.test(
  datboth$LofStd[datboth$Category == "Cahi"],
  datboth$LofStd[datboth$Category == "Caib"]
)
cat("\n=== Comparison of standardized LoF between domestic goat and ibex ===\n")
print(lof_test_result)

# ============================================================================
# Visualization 1: Total LoF variants
# ============================================================================

ggplot(datboth[datboth$Category %in% c("Cahi", "CahiNew", "hybrids", "Caib"), ]) + 
  geom_col(aes(x = reorder(V1, -LofTot), y = LofTot, fill = filt)) + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  facet_grid(. ~ Category, space = "free", scale = "free") + 
  ylim(0, 125) +
  labs(title = "Total LoF variants",
       x = "Individual ID",
       y = "Total LoF count",
       fill = "Reference")
ggsave("HighTot.stacked.pdf", width = 12, height = 6)

# ============================================================================
# Visualization 2: Heterozygous LoF variants (Figure 4D)
# ============================================================================

ggplot(datboth[datboth$Category %in% c("Cahi", "CahiNew", "hybrids", "Caib"), ]) + 
  geom_col(aes(x = reorder(V1, -LofHet), y = LofHet, fill = filt)) + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  facet_grid(. ~ Category, space = "free", scale = "free") +
  labs(title = "Heterozygous LoF variants (Figure 4D)",
       x = "Individual ID",
       y = "Heterozygous LoF count",
       fill = "Reference")
ggsave("Fig4D_LofHet.stacked.pdf", width = 12, height = 6)

# ============================================================================
# Visualization 3: Standardized LoF load - Boxplot comparison
# ============================================================================

ggplot(datboth[datboth$Category %in% c("Cahi", "Caib"), ]) + 
  geom_boxplot(aes(x = Category, y = LofTot / ModifTot, fill = filt)) + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  facet_grid(. ~ Category, space = "free", scale = "free") +
  labs(title = "Standardized LoF load (LofTot/ModifTot)",
       x = "Species category",
       y = "Standardized LoF load",
       fill = "Reference")
ggsave("LofTot.Std.box.onlyCahiib.pdf", width = 5, height = 5)

# ============================================================================
# Analysis 3: Combined analysis across both reference states
# ============================================================================
# Calculate combined LoF loads across both reference comparisons

datiibex$lofsBOTH <- datboth$LofTot[datboth$filt == "WOdomFixRef"] + 
                     datboth$LofTot[datboth$filt == "WOcaibFixRef"]
datiibex$modifsBOTH <- datboth$ModifTot[datboth$filt == "WOdomFixRef"] + 
                       datboth$ModifTot[datboth$filt == "WOcaibFixRef"]
datiibex$lofsHetBOTH <- datboth$LofHet[datboth$filt == "WOdomFixRef"] + 
                        datboth$LofHet[datboth$filt == "WOcaibFixRef"]

# ============================================================================
# Summary Statistics: LoF heterozygous variants across categories
# ============================================================================

cat("\n=== LoF Heterozygous Variants Summary Statistics ===\n")
cat("\nHybrids:\n")
print(summary(datiibex$lofsHetBOTH[datiibex$Category == "hybrids"]))

cat("\nIbex (Caib):\n")
print(summary(datiibex$lofsHetBOTH[datiibex$Category == "Caib"]))

cat("\nDomestic goat (Cahi & CahiNew):\n")
print(summary(datiibex$lofsHetBOTH[datiibex$Category %in% c("Cahi", "CahiNew")]))

# ============================================================================
# Visualization 4: Additional figures
# ============================================================================

# Stacked bar plot for all samples
ggplot(datboth[datboth$Category %in% c("Cahi", "CahiNew", "hybrids", "Caib"), ]) + 
  geom_col(aes(x = reorder(V1, -LofHet), y = LofHet, fill = filt)) + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  facet_grid(. ~ Category, space = "free", scale = "free") + 
  ylim(0, 125) +
  labs(title = "Heterozygous LoF variants by individual",
       x = "Individual ID",
       y = "Heterozygous LoF count",
       fill = "Reference")
ggsave("LofHet.stacked.pdf", width = 12, height = 6)

# Supplementary Figure: LoF standardized load in hybrids
ggplot(datiibex[datiibex$Category %in% c("hybrids"), ]) + 
  geom_col(aes(x = reorder(IndID, -PropGoat), y = lofsBOTH / modifsBOTH)) + 
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(title = "Standardized LoF load in hybrid individuals",
       x = "Individual ID (ordered by goat ancestry proportion)",
       y = "Standardized LoF load")
ggsave("Fig_Sx_LofStd.hybrids.pdf", width = 5, height = 4)

# Relationship between goat ancestry and LoF load in hybrids
ggplot(datiibex[datiibex$Category %in% c("hybrids"), ]) + 
  geom_point(aes(x = 1 - PropGoat, y = lofsBOTH / modifsBOTH), size = 2) + 
  theme_light() +
  labs(title = "LoF load vs. genetic ancestry in hybrids",
       x = "Proportion ibex ancestry (1 - PropGoat)",
       y = "Standardized LoF load")
ggsave("LofStd.PropGoat.points.pdf", width = 3, height = 4)

ggplot(datiibex[datiibex$Category %in% c("Cahi", "CahiNew", "hybrids","Caib"),]) + geom_col(aes(x=reorder(V1,PropGoat),y=lofsBOTH/modifsBOTH, fill=filt)) + facet_grid(.~Category, scale="free", space = "free")

#***** Figure 4C *********
ggplot(datiibex[datiibex$Category %in% c("Cahi", "CahiNew", "hybrids","Caib"),]) + geom_boxplot(aes(x=Category,y=lofsBOTH/modifsBOTH, fill=Category)) + theme_light()
ggsave("Fig4C_LofStd.allthree.box.pdf", width=5, height=5)
#*****************

##################################################
# ***  heterosis canditates
##################################################
dat <- read.table(file.path(DATA_DIR, "ibex.HIGH.CaibFixDerived.EXPR.het.scount"))
dat$HighHet <- dat$V2
# Note: homref comparison is interpreted as derived alleles fixed in ibex population
dat$HighHomDerived <- read.table(file.path(DATA_DIR, "ibex.HIGH.CaibFixDerived.EXPR.homref.scount"))$V2

dat$LofHet <- read.table(file.path(DATA_DIR, "ibex.LOF.CaibFixDerived.EXPR.het.scount"))$V2
dat$LofHomDerived <- read.table(file.path(DATA_DIR, "ibex.LOF.CaibFixDerived.EXPR.homref.scount"))$V2

dat$ModifHet <- read.table(file.path(DATA_DIR, "ibex.MODIFIER.CaibFixDerived.EXPR.het.scount"))$V2
dat$ModifHomDerived <- read.table(file.path(DATA_DIR, "ibex.MODIFIER.CaibFixDerived.EXPR.homref.scount"))$V2

dat$HighTot <- dat$HighHet + 2 * dat$HighHomDerived
dat$LofTot <- dat$LofHet + 2 * dat$LofHomDerived
dat$ModifTot <- dat$ModifHet + 2 * dat$ModifHomDerived

dati <- merge(dat, info, by.x = "V1", by.y = "Accession")
dati$Category <- factor(dati$Category, 
                        levels = c("Cahi", "hybrids", "Caib", "Capy", "other"))

# ============================================================================
# Visualization 7: LoF in heterosis candidates
# ============================================================================

ggplot(dati[dati$Category %in% c("hybrids"), ]) + 
  geom_col(aes(x = reorder(IndID, -PropGoat), y = LofTot, fill = "blue")) + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  facet_grid(. ~ Category, space = "free", scale = "free") +
  labs(title = "Total LoF variants in heterosis candidates",
       x = "Individual ID (ordered by goat ancestry)",
       y = "Total LoF count")
ggsave("LofTot.heterosis.stacked.pdf", width = 3, height = 3)

# ============================================================================
# Visualization 8: Heterozygous LoF in heterosis candidates - Figure 4B
# ============================================================================

ggplot(dati[dati$Category %in% c("hybrids"), ]) + 
  geom_col(aes(x = reorder(IndID, -PropGoat), y = LofHet, fill = "blue")) + 
  theme_light() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  facet_grid(. ~ Category, space = "free", scale = "free") +
  labs(title = "Heterozygous LoF variants in heterosis candidates (Figure 4B)",
       x = "Individual ID (ordered by goat ancestry)",
       y = "Heterozygous LoF count")
ggsave("Fig4B_LofHet.heterosis.stacked.pdf", width = 3, height = 3)

# ============================================================================
# END OF ANALYSIS
# ============================================================================

