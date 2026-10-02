
#' Separation index: predict gene-set null inflation before testing
#'
#' Computes the median absolute per-gene t-statistic for a grouping and
#' predicts the standardized empirical null floor from it. Derived across 17
#' cohort-grouping combinations in TCGA (Spearman 0.902 across 25 combinations, floor = 2.37 x
#' median|t|^1.38, log-log R2 = 0.83). Random groupings give median|t| near 0.6 and floors
#' near the theoretical 1.96; strong tumour-state groupings give median|t|
#' above 2 and floors above 5.
#'
#' @param expr Numeric matrix, genes in rows, samples in columns.
#' @param group Grouping variable of length ncol(expr).
#' @return A one-row data frame: median_abs_t, predicted_floor,
#'   ratio_to_theory and a verdict.
#' @export
sep_index <- function(expr, group) {
  if (!requireNamespace("limma", quietly = TRUE))
    stop("sep_index requires the limma package")
  expr <- expr[apply(expr, 1, stats::sd) > 0, , drop = FALSE]
  fit <- limma::eBayes(limma::lmFit(expr, stats::model.matrix(~ group)))
  tt <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")$t
  mt <- stats::median(abs(tt))
  data.frame(median_abs_t = round(mt, 3),
             predicted_floor = round(2.37 * mt^1.38, 2),
             ratio_to_theory = round(2.37 * mt^1.38 / 1.96, 2),
             verdict = if (mt < 0.9) "calibrated"
                       else if (mt < 1.7) "moderate inflation"
                       else "severe inflation",
             stringsAsFactors = FALSE)
}
