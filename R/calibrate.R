#' Mean inter-gene correlation of a gene set
#'
#' @param expr Numeric matrix, genes in rows (named), samples in columns.
#' @param genes Character vector of gene identifiers.
#' @param group Optional grouping vector. When supplied, correlation is computed
#'   on residuals after removing the group effect, so a real difference between
#'   groups does not inflate the estimate.
#' @return Mean pairwise correlation, or NA if fewer than three genes are present.
#' @export
set_coherence <- function(expr, genes, group = NULL) {
  g <- intersect(genes, rownames(expr))
  if (length(g) < 3) return(NA_real_)
  M <- expr[g, , drop = FALSE]
  if (!is.null(group)) M <- t(stats::resid(stats::lm(t(M) ~ group)))
  r <- stats::cor(t(M))
  mean(r[upper.tri(r)])
}

#' Empirical calibration for per-sample gene set scores
#'
#' Compares a gene set's association with an outcome against random gene sets
#' matched on size and mean expression decile, and reports how far the resulting
#' null sits from the parametric one.
#'
#' @section What this tests:
#' The empirical p-value answers a competitive question: is this set more
#' associated with the outcome than a typical gene set of the same size and
#' expression level? The nominal p-value answers the self-contained question of
#' whether the set is associated at all (Goeman and Buhlmann 2007). Where
#' background differential expression is high the two diverge sharply. In TCGA
#' lung adenocarcinoma, 55 percent of genes differ individually between
#' LKB1-deficient and intact tumours; 50 of 50 Hallmark sets are self-contained
#' significant and 9 are competitively significant.
#'
#' @section Known limitation:
#' This null is computed within one fixed grouping, so it captures variation
#' from which genes were drawn and not from which samples fall in each group.
#' For a real gene set the second component is larger: on HALLMARK_E2F_TARGETS
#' in TCGA LUAD the within-split null has SD 0.047 while the set's own
#' coefficient varies across random splits with SD 0.110. Type I error against
#' random sample splits reaches 0.30 for some sets.
#'
#' Two repairs were tested and neither works. Matching draws on inter-gene
#' correlation is impossible: real programs reach mean correlation 0.245 while
#' random draws of the same size reach at most 0.021. Combining gene
#' randomization with sample permutation reverts to the self-contained null
#' (Maciejewski 2014).
#'
#' Treat p_empirical as a diagnostic, not a corrected p-value. For inference
#' under dependence use limma::camera, limma::roast, or rSEA::SEA.
#'
#' @param expr Numeric matrix of log-scale expression, genes in rows (rownames
#'   required), samples in columns.
#' @param sets Named list of character vectors.
#' @param group Grouping variable of length ncol(expr) for model = "lm".
#'   Ignored when model = "cox".
#' @param covariates Optional data frame of covariates, one row per sample.
#' @param model "lm" for a linear model on group, or "cox" for Cox
#'   proportional hazards on time and event.
#' @param time Survival time, required when model = "cox".
#' @param event Event indicator, required when model = "cox".
#' @param n_null Random sets per tested set. Default 100; p-value resolution is
#'   1/(n_null + 1).
#' @param score_fn Function of (expr, sets) returning one row per set. Defaults
#'   to the mean of standardized expression.
#' @param seed Random seed.
#'
#' @return A data frame with one row per set: beta, p_nominal, p_empirical,
#'   floor_95 (95th percentile of the absolute null coefficient), floor_std
#'   (that floor in units of its own standard error, comparable to 1.96) and
#'   rho_set (residual coherence).
#'
#' @references
#' Goeman JJ, Buhlmann P (2007). Analyzing gene expression data in terms of gene
#' sets: methodological issues. Bioinformatics 23, 980-987.
#'
#' Wu D, Smyth GK (2012). Camera: a competitive gene set test accounting for
#' inter-gene correlation. Nucleic Acids Research 40, e133.
#'
#' Maciejewski H (2014). Gene set analysis methods: statistical models and
#' methodological differences. Briefings in Bioinformatics 15, 504-518.
#'
#' Phipson B, Smyth GK (2010). Permutation P-values should never be zero.
#' Statistical Applications in Genetics and Molecular Biology 9, 39.
#'
#' @examples
#' set.seed(1)
#' e <- matrix(rnorm(2000 * 40), 2000, 40,
#'             dimnames = list(paste0('g', 1:2000), paste0('s', 1:40)))
#' grp <- rep(c(TRUE, FALSE), each = 20)
#' calibrate_geneset(e, list(setA = paste0('g', 1:50)), grp, n_null = 50)
#' @export
calibrate_geneset <- function(expr, sets, group = NULL, covariates = NULL,
                              model = c("lm", "cox"), time = NULL, event = NULL,
                              n_null = 100, score_fn = NULL, seed = 1) {
  model <- match.arg(model)
  stopifnot(is.matrix(expr), !is.null(rownames(expr)))
  if (is.null(names(sets))) stop("sets must be a named list")

  if (model == "lm") {
    if (is.null(group)) stop("group is required when model = 'lm'")
    if (length(group) != ncol(expr)) stop("group must have length ncol(expr)")
  } else {
    if (!requireNamespace("survival", quietly = TRUE))
      stop("model = 'cox' requires the survival package")
    if (is.null(time) || is.null(event))
      stop("time and event are required when model = 'cox'")
    if (length(time) != ncol(expr) || length(event) != ncol(expr))
      stop("time and event must have length ncol(expr)")
  }
  if (!is.null(covariates) && nrow(covariates) != ncol(expr))
    stop("covariates must have one row per sample")

  expr <- expr[apply(expr, 1, stats::sd) > 0, , drop = FALSE]
  if (nrow(expr) < 100) stop("too few variable genes in expr")

  if (is.null(score_fn)) {
    zz <- t(scale(t(expr)))
    score_fn <- function(m, s) t(vapply(s, function(g) {
      g <- intersect(g, rownames(zz))
      colMeans(zz[g, , drop = FALSE])
    }, numeric(ncol(zz))))
  }

  zsc <- function(x) (x - mean(x)) / stats::sd(x)
  gm  <- rowMeans(expr)
  dec <- cut(gm, stats::quantile(gm, 0:10 / 10), include.lowest = TRUE, labels = FALSE)
  names(dec) <- rownames(expr)
  univ <- setdiff(rownames(expr), unlist(sets))
  pool <- split(univ, dec[univ])

  cvn <- if (is.null(covariates)) character(0) else names(covariates)
  base <- if (is.null(covariates)) data.frame(row.names = seq_len(ncol(expr))) else covariates

  fitb <- function(y) {
    d <- base; d$.y <- zsc(y)
    if (model == "lm") {
      d$.grp <- group
      f <- stats::as.formula(paste(".y ~", paste(c(".grp", cvn), collapse = " + ")))
      stats::coef(summary(stats::lm(f, data = d)))[2, c(1, 4)]
    } else {
      d$.t <- time; d$.e <- event
      f <- stats::as.formula(paste("survival::Surv(.t, .e) ~",
                                   paste(c(".y", cvn), collapse = " + ")))
      s <- summary(survival::coxph(f, data = d))$coefficients
      s[".y", c(1, 5)]
    }
  }

  se_fac <- if (model == "lm") {
    sqrt(1 / sum(group == group[1]) + 1 / sum(group != group[1]))
  } else 1 / sqrt(sum(event == 1))

  set.seed(seed)
  out <- lapply(names(sets), function(nm) {
    g <- intersect(sets[[nm]], rownames(expr))
    if (length(g) < 5) return(NULL)
    obs <- fitb(score_fn(expr, list(.obs = g))[1, ])

    need  <- table(dec[g])
    draws <- replicate(n_null, unlist(lapply(names(need), function(k)
      sample(pool[[k]], min(need[[k]], length(pool[[k]]))))), simplify = FALSE)
    S  <- score_fn(expr, stats::setNames(draws, paste0("n", seq_len(n_null))))
    nb <- vapply(seq_len(nrow(S)), function(i) unname(fitb(S[i, ])[1]), numeric(1))

    data.frame(
      set = nm, m = length(g), model = model,
      beta = unname(obs[1]), p_nominal = unname(obs[2]),
      p_empirical = (sum(abs(nb) >= abs(obs[1])) + 1) / (n_null + 1),
      floor_95 = unname(stats::quantile(abs(nb), 0.95)),
      floor_std = unname(stats::quantile(abs(nb), 0.95)) / se_fac,
      rho_set = set_coherence(expr, g, if (model == "lm") group else NULL),
      stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}
