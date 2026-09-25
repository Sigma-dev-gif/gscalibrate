#' Mean inter-gene correlation of a gene set
#'
#' @param expr Numeric matrix, genes in rows (named), samples in columns.
#' @param genes Character vector of gene identifiers.
#' @return Mean absolute-signed pairwise correlation, or NA if fewer than three
#'   genes are present in \code{expr}.
#' @export
set_coherence <- function(expr, genes) {
  g <- intersect(genes, rownames(expr))
  if (length(g) < 3) return(NA_real_)
  r <- stats::cor(t(expr[g, , drop = FALSE]))
  mean(r[upper.tri(r)])
}

#' Empirical calibration for per-sample gene set scores
#'
#' Computes an empirical p-value for the association between a gene set score
#' and a grouping variable, using a null distribution built from random gene
#' sets matched on set size and mean expression decile.
#'
#' @section Important limitation:
#' Draw-based nulls are anti-conservative when the tested set is internally
#' correlated, because random draws cannot reproduce the coherence of a real
#' gene set. In TCGA lung adenocarcinoma, tissue programs reach mean inter-gene
#' correlation of 0.245 while random draws of equal size reach at most 0.021.
#' Simulated type I error for a coherent set reaches 0.125 at 250 genes and
#' 0.255 at 900 genes. This function reports the coherence gap and warns when
#' it is large. For coherent sets, prefer \code{limma::roast} or
#' \code{limma::cameraPR}, which estimate correlation from the set itself.
#'
#' @param expr Numeric matrix of expression on a log scale, genes in rows
#'   (rownames required), samples in columns.
#' @param sets Named list of character vectors.
#' @param group Grouping variable, length \code{ncol(expr)}.
#' @param covariates Optional data frame of covariates, one row per sample.
#' @param n_null Number of random sets per tested set. Default 100; resolution
#'   of the empirical p-value is \code{1/n_null}.
#' @param score_fn Function of \code{(expr, sets)} returning a score matrix with
#'   one row per set. Defaults to mean of standardized expression.
#' @param seed Random seed.
#' @param coherence_warn Warn when the tested set's coherence exceeds the
#'   maximum achieved by random draws times this factor. Default 2.
#'
#' @return A data frame with one row per set: \code{beta}, \code{p_nominal},
#'   \code{p_empirical}, \code{floor_95} (95th percentile of the absolute null
#'   coefficient), \code{floor_std} (that floor in units of its own standard
#'   error; compare to 1.96), \code{rho_set}, \code{rho_null_max}, and
#'   \code{reliable} (FALSE when the coherence gap is large).
#'
#' @references
#' Wu D, Smyth GK (2012). Camera: a competitive gene set test accounting for
#' inter-gene correlation. \emph{Nucleic Acids Research} 40, e133.
#'
#' @examples
#' set.seed(1)
#' e <- matrix(rnorm(2000 * 40), 2000, 40,
#'             dimnames = list(paste0('g', 1:2000), paste0('s', 1:40)))
#' grp <- rep(c(TRUE, FALSE), each = 20)
#' calibrate_geneset(e, list(setA = paste0('g', 1:50)), grp, n_null = 50)
#' @export
calibrate_geneset <- function(expr, sets, group, covariates = NULL,
                              n_null = 100, score_fn = NULL, seed = 1,
                              coherence_warn = 2) {
  stopifnot(is.matrix(expr), !is.null(rownames(expr)))
  if (length(group) != ncol(expr))
    stop('group must have length ncol(expr)')
  if (!is.null(covariates) && nrow(covariates) != ncol(expr))
    stop('covariates must have one row per sample')
  if (is.null(names(sets))) stop('sets must be a named list')

  expr <- expr[apply(expr, 1, stats::sd) > 0, , drop = FALSE]
  if (nrow(expr) < 100) stop('too few variable genes in expr')

  if (is.null(score_fn)) {
    zz <- t(scale(t(expr)))
    score_fn <- function(m, s) t(vapply(s, function(g) {
      g <- intersect(g, rownames(zz))
      colMeans(zz[g, , drop = FALSE])
    }, numeric(ncol(zz))))
  }

  zsc <- function(x) (x - mean(x)) / stats::sd(x)
  resid_coh <- function(g) {
    if (length(g) < 3) return(NA_real_)
    R <- t(stats::resid(stats::lm(t(expr[g, , drop = FALSE]) ~ group)))
    r <- stats::cor(t(R)); mean(r[upper.tri(r)])
  }
  gm  <- rowMeans(expr)
  dec <- cut(gm, stats::quantile(gm, 0:10 / 10),
                    include.lowest = TRUE, labels = FALSE)
  names(dec) <- rownames(expr)
  univ <- setdiff(rownames(expr), unlist(sets))
  pool <- split(univ, dec[univ])

  dat <- data.frame(.grp = group)
  if (!is.null(covariates)) dat <- cbind(dat, covariates)
  rhs <- paste(names(dat), collapse = ' + ')
  form <- stats::as.formula(paste('.y ~', rhs))

  set.seed(seed)
  out <- lapply(names(sets), function(nm) {
    g <- intersect(sets[[nm]], rownames(expr))
    if (length(g) < 5) return(NULL)
    need <- table(dec[g])
    draws <- replicate(n_null, unlist(lapply(names(need), function(k)
      sample(pool[[k]], min(need[[k]], length(pool[[k]]))))), simplify = FALSE)

    S    <- score_fn(expr, c(list(.obs = g), stats::setNames(draws, paste0('n', seq_len(n_null)))))
    fitb <- function(y) {
      d <- dat; d$.y <- zsc(y)
      stats::coef(summary(stats::lm(form, data = d)))[2, c(1, 4)]
    }
    obs <- fitb(S[1, ])
    nb  <- vapply(2:nrow(S), function(i) fitb(S[i, ])[1], numeric(1))

    se_fac <- sqrt(1 / sum(group == group[1]) + 1 / sum(group != group[1]))
    rho_s  <- resid_coh(g)
    rho_n  <- max(vapply(draws, function(d) set_coherence(expr, d), numeric(1)),
                  na.rm = TRUE)

    data.frame(set = nm, m = length(g),
               beta = unname(obs[1]), p_nominal = unname(obs[2]),
               p_empirical = (sum(abs(nb) >= abs(obs[1])) + 1) / (length(nb) + 1),
               floor_95 = unname(stats::quantile(abs(nb), 0.95)),
               floor_std = unname(stats::quantile(abs(nb), 0.95)) / se_fac,
               rho_set = rho_s, rho_null_max = rho_n,
               recommended = if (length(g) >= 150) 'roast' else 'empirical',
               reliable = !is.na(rho_s) && rho_s <= coherence_warn * rho_n,
               stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, out)
  if (any(!res$reliable))
    warning('Sets with reliable = FALSE are more internally correlated than any ',
            'random draw achieved; the empirical null is anti-conservative for ',
            'these. Use limma::roast or limma::cameraPR instead.', call. = FALSE)
  rownames(res) <- NULL
  res
}
