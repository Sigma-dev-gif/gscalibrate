## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## -----------------------------------------------------------------------------
library(gscalibrate)
set.seed(1)

n_genes <- 3000; n_samples <- 80
expr <- matrix(rnorm(n_genes * n_samples), n_genes, n_samples,
               dimnames = list(paste0("gene", seq_len(n_genes)),
                               paste0("s", seq_len(n_samples))))
group <- rep(c(TRUE, FALSE), each = n_samples / 2)

sets <- list(setA = paste0("gene", 1:120),
             setB = paste0("gene", 201:320))

calibrate_geneset(expr, sets, group, n_null = 100)

## -----------------------------------------------------------------------------
set.seed(2)
expr2 <- expr

# 40% of all genes differ slightly between groups
de <- sample(n_genes, 0.4 * n_genes)
expr2[de, group] <- expr2[de, group] + rnorm(length(de), 0, 0.3)

# setA shares a latent factor: a coherent set
lat <- rnorm(n_samples)
idx <- match(sets$setA, rownames(expr2))
expr2[idx, ] <- sqrt(0.15) * matrix(rep(lat, each = length(idx)),
                                    length(idx), n_samples) +
                sqrt(0.85) * expr2[idx, ]

res <- calibrate_geneset(expr2, sets, group, n_null = 100)

        

## ----eval = FALSE-------------------------------------------------------------
# library(limma)
# design <- model.matrix(~ group)
# roast(expr2, index = which(rownames(expr2) %in% sets$setA),
#       design = design, contrast = 2)
# cameraPR(statistic, index)   # direct replacement for pre-ranked GSEA

