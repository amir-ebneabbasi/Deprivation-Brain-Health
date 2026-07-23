# Deprivation–Brain–Disease Mediation Analysis

## Overview

This repository contains the code used to quantify the neuroanatomical mediation of neighbourhood deprivation on psychiatric and neurological disease risk.

The analysis was developed for two large population-based cohorts:

* **UK Biobank (UKB):** approximately 500,000 adults aged 40–69 years
* **Adolescent Brain Cognitive Development (ABCD) Study:** 11,878 children aged 9–10 years

Across both cohorts, mediation analyses were performed between neighbourhood deprivation, regional brain imaging phenotypes, and psychiatric and neurological disorders.

---

## Statistical framework

For each deprivation–brain–disease combination, two regression models are fitted.

### 1. Mediator model

Linear regression:

```
Brain volume ~ Deprivation + Covariates
```

The coefficient for deprivation represents **path a**.

### 2. Outcome model

Binomial logistic regression:

```
Disease ~ Deprivation + Brain volume + Covariates
```

The coefficient for brain volume represents **path b**, while the coefficient for deprivation represents the **direct effect**.

The indirect (mediation) effect is calculated as:

```
Indirect effect = a × b
```

Bootstrap resampling is used to estimate:

* indirect effects
* percentile confidence intervals
* empirical two-sided p-values

---

## Features

* Bootstrap mediation analysis
* Linear regression for continuous brain imaging phenotypes
* Binomial logistic regression for binary disease outcomes
* Automatic control-group construction
* Reproducible bootstrap sampling
* Parallel processing of thousands of mediation models
* SLURM array support for high-performance computing

---

## License

This repository is released under the MIT License.
