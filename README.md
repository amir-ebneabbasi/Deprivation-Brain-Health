# Deprivation–Brain–Disease and Genetic Analyses

## Overview

This repository contains code for analysing relationships between neighbourhood deprivation, brain phenotypes, disease risk, and genetic variation in large neuroimaging cohorts, including UK Biobank (UKB) and the Adolescent Brain Cognitive Development (ABCD) Study.

The repository includes two main components:

1. **Deprivation–brain–disease mediation analysis**
2. **Genetic analysis using KING and GENESIS**

---

## Deprivation–Brain–Disease Mediation

`dep_main.py` tests whether regional brain phenotypes mediate associations between neighbourhood deprivation and psychiatric or neurological disease.

Two models are fitted:

```text
Brain phenotype ~ Deprivation + Covariates
Disease ~ Deprivation + Brain phenotype + Covariates
```

The mediation effect is calculated as:

```text
Indirect effect = a × b
```

Bootstrap resampling is used to estimate indirect effects, confidence intervals, and empirical p-values.

---

## Genetic Analysis

The genetic pipeline performs relatedness estimation, population structure correction, genome-wide association testing, and heritability estimation.

### `plink_to_gds.R`

Converts merged autosomal PLINK files (`BED/BIM/FAM`) to GDS format using `SNPRelate`.

**Output:** `genotype_autosomes.gds`

### `king.sh`

Runs KING to estimate pairwise genetic relatedness up to third-degree relatives.

**Output:** `eur_king.kin0`

### `king_to_matrix.R`

Converts KING relatedness estimates into a kinship matrix compatible with `GENESIS`.

**Output:** `eur_kinship_matrix.rds`

### `PC_air.R`

Runs PC-AiR using LD-pruned SNPs and KING relatedness estimates to obtain ancestry principal components while accounting for related individuals.

**Output:** `mypc.rds`

### `PC_relate.R`

Runs PC-Relate using ancestry PCs and unrelated reference samples to estimate ancestry-adjusted genetic relatedness and construct a sparse genetic relationship matrix (GRM).

**Outputs:** `mypcrelate.rds` and `mypcrelate_sparse.rds`

### `genesis.R`

Runs the final GENESIS analysis for each brain phenotype using a SLURM array. It:

* fits a linear mixed model using the PC-Relate GRM
* adjusts for demographic, imaging, and PC-AiR covariates
* performs single-SNP genome-wide association testing
* estimates SNP-based heritability
* exports GWAS summary statistics

**Outputs:** `assoc_<ID>.rds`, `heritability_<ID>.rds`, and `GWAS_<ID>.txt`

---

## Genetic Workflow

```text
PLINK BED/BIM/FAM
       ↓
plink_to_gds.R
       ↓
KING → king_to_matrix.R
       ↓
PC_air.R
       ↓
PC_relate.R
       ↓
Sparse GRM
       ↓
genesis.R
       ↓
GWAS + Heritability
```

## License

This repository is released under the MIT License.
