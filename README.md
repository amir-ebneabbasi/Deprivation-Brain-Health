

# Deprivation–Brain–Disease and Genetic Analyses

**Preprint:** Ebneabbasi A, Warrier V, Montagnese M, Romero Garcia R, Bethlehem RAI, Rittman T. *Mapping the Health Burden of Neighbourhood Deprivation: Neurobiological Evidence Across the Life Span.* medRxiv (2026). [https://doi.org/10.64898/2026.08.29.26361714](https://doi.org/10.64898/2026.08.29.26361714) · [Preprint page](https://www.medrxiv.org/content/10.64898/2026.08.29.26361714v1) · [PDF](https://www.medrxiv.org/content/10.64898/2026.08.29.26361714v1.full.pdf)

## Overview

This repository contains code for analysing relationships between neighbourhood deprivation, brain phenotypes, disease risk across three cohorts spanning the life span:

| Cohort | Sample | Age range |
|---|---|---|
| HEALthy Brain and Child Development (HBCD) Study | n = 84 | 0–4 weeks postnatal |
| Adolescent Brain Cognitive Development (ABCD) Study | n = 4,792 | 9–10 years |
| UK Biobank (UKB) | ~500,000 adults | 44–87 years |

The repository includes two main components:

1. **Deprivation–brain–disease mediation analysis** (`dep_main.py`)
2. **Genetic analysis using KING and GENESIS** (R and shell scripts)

---

## Table of contents

- [Computing environment](#computing-environment)
- [Software dependencies](#software-dependencies)
- [Data availability](#data-availability)
- [Data processing and provenance](#data-processing-and-provenance)
- [Deprivation–brain–disease mediation](#deprivationbraindisease-mediation)
- [Genetic analysis](#genetic-analysis)
- [Running on SLURM](#running-on-slurm)
- [Citation](#citation)
- [License](#license)

---

## Computing environment

All analyses were run on a high-performance computing (HPC) cluster.

---

## Software versions

### Versions used in the study

| Software | Version | Used for |
|---|---|---|
| Python | 3.11 | Mediation and statistical analyses |
| statsmodels | 0.14.4 | Regression models |
| SNPRelate | 1.34.1 | GDS conversion and LD pruning |
| KING | 2.3.2 | Pairwise relatedness (up to third degree) |
| GENESIS | 2.30.0 | Ancestry PCs, kinship, linear mixed models |
| PC-AiR | 0.8.0 | Ancestry principal components accounting for relatedness |
| PC-Relate | 1.0.0 | Ancestry-adjusted GRM |

PC-AiR and PC-Relate are run through the GENESIS package (`pcair()` and `pcrelate()`).

Multiple-comparison correction used the Benjamini–Hochberg false discovery rate (FDR) procedure.


### Example environment setup

```bash
pip install statsmodels==0.14.4 numpy pandas scipy

# R (genetic analysis) — run inside R
# if (!require("BiocManager")) install.packages("BiocManager")
# BiocManager::install(c("SNPRelate", "GENESIS"))   # target GENESIS 2.30.0, SNPRelate 1.34.1
```

KING (v2.3.2) is a standalone binary; download it from the [KING website](https://www.kingrelatedness.com/) and make sure it is on your `PATH` (or loaded as a module on your cluster).

### Installation time

Setting up the environment is quick. On a standard HPC node with internet access, installing the Python packages takes about 1–2 minutes, downloading the KING binary takes under a minute, and installing the R/Bioconductor packages (GENESIS and SNPRelate, plus their dependencies) takes roughly 5–10 minutes.

---

## Data availability

Participant-level data are controlled-access and cannot be redistributed by the authors. This repository therefore contains code only.

- **ABCD and HBCD:** available to eligible researchers through the [NIH Brain Development Cohorts Data Hub](https://www.nbdc-datahub.org/), subject to approval of a Data Use Certification and completion of the required training.
- **UK Biobank:** available to eligible researchers through the [UK Biobank access process](https://www.ukbiobank.ac.uk/use-our-data/apply-for-access/).

---

## Data processing and provenance

- **HBCD:** imaging and genetic data were processed by the HBCD Study, and preprocessed data were downloaded for the present analyses.
- **UKB genetic data:** processed by the UK Biobank team.
- **ABCD and UKB imaging:** processed using FreeSurfer (v6.0.1) workflows, as described at:
  - ABCD: [ucam-department-of-psychiatry/ABCD](https://github.com/ucam-department-of-psychiatry/ABCD)
  - UKB: [ucam-department-of-psychiatry/UKB](https://github.com/ucam-department-of-psychiatry/UKB)
- **ABCD genetic quality control:** followed the pipeline at [vwarrier/ABCD_geneticQC](https://github.com/vwarrier/ABCD_geneticQC).
- **Spin permutation analyses:** performed with code at [amir-ebneabbasi/Spatial-Colocation](https://github.com/amir-ebneabbasi/Spatial-Colocation).

---

## Deprivation–brain–disease mediation

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

where *a* is the deprivation coefficient in the brain model and *b* is the brain-phenotype coefficient in the disease model.

Bootstrap resampling is used to estimate indirect effects, confidence intervals, and empirical p-values. Multiple comparisons are controlled with Benjamini–Hochberg FDR.

---

## Genetic analysis

```text
plink_to_gds.R → ld_pruning.R → king.sh → king_to_matrix.R → pc_air.R → pc_relate.R → genesis.R
```

### `plink_to_gds.R`

Converts merged autosomal PLINK files (`BED/BIM/FAM`) to GDS format using `SNPRelate`.

**Output:** `genotype_autosomes.gds`

### `ld_pruning.R`

Performs linkage disequilibrium (LD) pruning of the autosomal genotype data using `SNPRelate`. Correlation-based pruning is applied to reduce redundancy among highly correlated genetic variants.

The resulting independent SNP set is used for downstream PC-AiR and PC-Relate analyses.

**Output:** `pruned_snps.txt`

### `king.sh`

Runs KING (v2.3.2) to estimate pairwise genetic relatedness up to third-degree relatives.

**Output:** `eur_king.kin0`

### `king_to_matrix.R`

Converts KING relatedness estimates into a kinship matrix compatible with `GENESIS`.

**Output:** `eur_kinship_matrix.rds`

### `pc_air.R`

Runs PC-AiR using LD-pruned SNPs and KING relatedness estimates to obtain ancestry principal components while accounting for related individuals.

**Output:** `pc.rds`

### `pc_relate.R`

Runs PC-Relate using ancestry PCs and unrelated reference samples to estimate ancestry-adjusted genetic relatedness and construct a sparse genetic relationship matrix (GRM).

**Outputs:** `pcrelate.rds` and `pcrelate_sparse.rds`

### `genesis.R`

Runs the final GENESIS analysis for each brain phenotype using a SLURM array. It:

- fits a linear mixed model using the PC-Relate GRM
- adjusts for demographic, imaging, and PC-AiR covariates

**Outputs:** `kinship-adjusted estimates`

---

## Running on SLURM

TBC

---

## Citation

If you use this code, please cite:

```bibtex
@article{ebneabbasi2026deprivation,
  title   = {Mapping the Health Burden of Neighbourhood Deprivation: Neurobiological Evidence Across the Life Span},
  author  = {Ebneabbasi, Amir and Warrier, Varun and Montagnese, Marcella and Romero Garcia, Rafael and Bethlehem, Richard A.I. and Rittman, Timothy},
  journal = {medRxiv},
  year    = {2026},
  doi     = {10.64898/2026.08.29.26361714}
}
```

---

## License

This repository is released under the MIT License.
