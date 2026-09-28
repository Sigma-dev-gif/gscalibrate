# gscalibrate

Empirical calibration for per-sample gene set scores.

## The problem

Per-sample gene set scores (ssGSEA, GSVA, mean-z) are routinely entered as
variables in regression and survival models with no correction for inter-gene
correlation. In a pre-registered pan-cancer analysis, 28 of 43 associations
significant after Benjamini-Hochberg correction (65%) failed an empirical null
built from matched-random gene sets.

## Install

```r
remotes::install_github('Sigma-dev-gif/gscalibrate')
```

## Use

```r
res <- calibrate_geneset(expr, sets, group, covariates = cov)
res[, c('set','beta','p_nominal','p_uncorrected','p_empirical','ratio')]
```

`p_uncorrected` is the naive empirical p-value; `p_empirical` is corrected.
`ratio` is how much too narrow the uncorrected null was — values above about 1.5
mean the naive p-value cannot be trusted.

## What the correction does

A gene-randomization null is computed within a fixed sample split, so it captures
only the variation from which genes were drawn. A real gene set's statistic also
varies with which samples land in each group, and that component is larger. In
TCGA lung adenocarcinoma, HALLMARK_E2F_TARGETS has a within-split null SD of 0.047
while its own coefficient varies across random splits with SD 0.110 — the null is
2.3x too narrow. Uncorrected type I error reaches 0.30.

`calibrate_geneset()` estimates the across-split spread directly, then centres and
rescales the null to match. Type I error across three cohorts and several gene sets:

| | uncorrected | corrected |
|---|---|---|
| COADREAD G2M checkpoint | 0.30 | 0.045 |
| BRCA G2M checkpoint | 0.19 | 0.055 |
| LUAD E2F targets | 0.15 | 0.043 |
| Notch signalling (already calibrated) | 0.04 | 0.055 |

The corrected p-values agree closely with `limma::roast`, which reaches the same
place by rotating residuals.

## Reference

Wu D, Smyth GK (2012). Camera: a competitive gene set test accounting for
inter-gene correlation. *Nucleic Acids Research* 40, e133.

Goeman JJ, Buhlmann P (2007). Analyzing gene expression data in terms of gene
sets: methodological issues. *Bioinformatics* 23, 980-987.
