# Deprivation–Brain–Disease Mediation Analysis

## Overview

This repository contains the code used to quantify the neuroanatomical mediation of neighbourhood deprivation on psychiatric and neurological disease risk.

The pipeline was developed for analyses in:

* **UK Biobank (UKB)** (~500,000 adults)
* **Adolescent Brain Cognitive Development (ABCD) Study** (11,878 children)

Across both cohorts, mediation analyses were performed between deprivation measures, regional brain imaging phenotypes, and multiple psychiatric and neurological disorders.

## Statistical framework

The mediation framework consists of:

1. Linear regression

[
\text{Brain Volume} = \beta_0 + a(\text{Deprivation}) + \mathbf{Covariates}
]

2. Logistic regression

[
\log\left(\frac{P(\text{Disease})}{1-P(\text{Disease})}\right)
==============================================================

\beta_0
+
c'(\text{Deprivation})
+
b(\text{Brain Volume})
+
\mathbf{Covariates}
]

The indirect (mediation) effect is computed as

[
a \times b
]

Bootstrap sampling is used to estimate:

* indirect effects
* percentile confidence intervals
* empirical two-sided p-values

## Features

* Logistic regression for binary disease outcomes
* Linear regression for imaging phenotypes
* Automatic control-group construction
* Reproducible bootstrap sampling
* Parallel processing of thousands of mediation models

## License

This repository is released under the MIT License.
