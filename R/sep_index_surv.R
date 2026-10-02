
#' Separation index for survival outcomes
#'
#' Cox version of \code{\link{sep_index}}. Computes the median absolute
#' per-gene Cox z-statistic and predicts the empirical null floor.
#' IMPORTANT: this returns median |z| only. The power law relating
#' separation to the null floor was calibrated on binary groupings and does
#' NOT transfer to survival outcomes: across three TCGA cohorts the observed
#' floors were 3.19, 1.76 and 3.10 against predictions of 2.97, 2.52 and 1.74,
#' with the ordering wrong. Use median |z| as a relative indicator only, and
#' estimate the floor empirically with calibrate_geneset(model = 'cox').
#'
#' @param expr Numeric matrix, genes in rows, samples in columns.
#' @param time Survival time.
#' @param event Event indicator.
#' @param n_genes Subsample this many genes for speed; NULL uses all.
#' @return A one-row data frame.
#' @export
sep_index_surv <- function(expr, time, event, n_genes = 4000) {
  if (!requireNamespace("survival", quietly = TRUE))
    stop("sep_index_surv requires the survival package")
  expr <- expr[apply(expr, 1, stats::sd) > 0, , drop = FALSE]
  k <- !is.na(time) & !is.na(event) & time > 0
  expr <- expr[, k, drop = FALSE]
  sv <- survival::Surv(time[k], event[k])
  g <- if (!is.null(n_genes) && n_genes < nrow(expr))
         sample(rownames(expr), n_genes) else rownames(expr)
  z <- vapply(g, function(i) {
    s <- try(summary(survival::coxph(sv ~ scale(expr[i, ])))$coefficients[1, 4],
             silent = TRUE)
    if (inherits(s, "try-error")) NA_real_ else s }, numeric(1))
  mt <- stats::median(abs(z), na.rm = TRUE)
  data.frame(median_abs_z = round(mt, 3),
             predicted_floor = round(2.37 * mt^1.38, 2),
             note = "floor prediction not calibrated for survival; see ?sep_index_surv",
             verdict = if (mt < 0.9) "low separation"
                       else if (mt < 1.7) "moderate separation"
                       else "high separation",
             stringsAsFactors = FALSE)
}
