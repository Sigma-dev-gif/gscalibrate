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
res[, c('set','beta','p_nominal','p_empirical','floor_std','reliable')]
```

`floor_std` is the 95th percentile of the null coefficient in units of its own
standard error. Compare it to 1.96: values well above that indicate the
parametric test is miscalibrated for this data and grouping.

## Important limitation

A draw-based null is **anti-conservative when the tested set is internally
coherent**, because random draws cannot reproduce the correlation structure that
defines a real gene set. Measured in TCGA lung adenocarcinoma: tissue programs
reach mean inter-gene correlation of 0.245, random draws of equal size at most
0.021. Simulated type I error reaches 0.125 at 250 genes and 0.255 at 900.

The `reliable` column flags sets where this applies. For those, use
`limma::roast` or `limma::cameraPR`, which estimate correlation from the set
itself rather than from draws.

## Reference

Wu D, Smyth GK (2012). Camera: a competitive gene set test accounting for
inter-gene correlation. *Nucleic Acids Research* 40, e133.
