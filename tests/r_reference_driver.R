#!/usr/bin/env Rscript
suppressMessages({
  library(condiments); library(jsonlite); library(slingshot)
})
args <- commandArgs(trailingOnly = TRUE)
fixture_path <- args[1]; output_path <- args[2]

set.seed(42)
toy <- condiments::create_differential_topology(n_cells = 500,
                                                shift = 10,
                                                unbalance_level = 0.9)
# toy is list(sd=DataFrame, mst=matrix). sd cols: Dim1 Dim2 lineages conditions.
df <- toy$sd
rd <- as.matrix(df[, c("Dim1", "Dim2")])
conds <- as.character(df$conditions)
lineage_init <- df$lineages   # in {-1, +1}; convert to {1, 2}
cluster_init <- ifelse(lineage_init == -1, "1", "2")

set.seed(42)
sds <- slingshot::slingshot(rd, clusterLabels = cluster_init, start.clus = "1")
pt <- slingshot::slingPseudotime(sds, na = FALSE)
cw <- slingshot::slingCurveWeights(sds)
cat("[ref] sds fitted; pt dims:", dim(pt), " cw dims:", dim(cw), "\n")

# Save sidecars for Python
write.csv(rd, sub("\\.rds$", "_reducedDim.csv", fixture_path), row.names = FALSE)
write.csv(pt, sub("\\.rds$", "_pseudotime.csv", fixture_path), row.names = FALSE)
write.csv(cw, sub("\\.rds$", "_cellWeights.csv", fixture_path), row.names = FALSE)
write.csv(data.frame(condition = conds),
          sub("\\.rds$", "_conditions.csv", fixture_path), row.names = FALSE)

# Run imbalance_score (smoothed)
set.seed(42)
imb <- condiments::imbalance_score(rd, conditions = conds, k = 10, smooth = 10)
cat("[ref] imbalance_score: length", length(imb$scaled_score), "\n")

# topologyTest
set.seed(42)
ttop <- tryCatch({
  condiments::topologyTest(sds = sds, conditions = conds, rep = 50)
}, error = function(e) { cat("[ref] topologyTest err:", e$message, "\n"); NULL })
top_pval <- if (!is.null(ttop)) min(as.numeric(ttop$p.value), na.rm = TRUE) else NA_real_

# progressionTest
set.seed(42)
ptest <- tryCatch({
  condiments::progressionTest(sds, conditions = conds)
}, error = function(e) { cat("[ref] progressionTest err:", e$message, "\n"); NULL })
prog_pval <- if (!is.null(ptest)) as.numeric(ptest$p.value)[1] else NA_real_
prog_stat <- if (!is.null(ptest)) as.numeric(ptest$statistic)[1] else NA_real_

# fateSelectionTest (= differentiationTest in older versions)
set.seed(42)
dtest <- tryCatch({
  condiments::fateSelectionTest(sds, conditions = conds)
}, error = function(e) { cat("[ref] fateSelectionTest err:", e$message, "\n"); NULL })
diff_pval <- if (!is.null(dtest)) as.numeric(dtest$p.value)[1] else NA_real_

out <- list(
  imbalance = list(
    score = as.numeric(imb$score),
    scaled_score = as.numeric(imb$scaled_score)
  ),
  topology = list(pvalue = top_pval),
  progression = list(pvalue = prog_pval, statistic = prog_stat),
  differentiation = list(pvalue = diff_pval)
)
write_json(out, output_path, auto_unbox = TRUE, digits = NA, na = "null", pretty = FALSE)
cat("[ref] wrote", output_path, "\n")
