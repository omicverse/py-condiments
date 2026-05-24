#!/usr/bin/env Rscript
# R reference plots for condiments vignette.
suppressPackageStartupMessages({
  library(condiments)
  library(slingshot)
  library(ggplot2)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Generate fixture
set.seed(42)
toy <- create_differential_topology(n_cells = 300, shift = 10)
df <- data.frame(toy$sd, conditions = toy$sd$conditions)

p <- ggplot(df, aes(x = Dim1, y = Dim2, col = conditions)) +
  geom_point(alpha = 0.5, size = 1.5) +
  theme_classic() +
  labs(x = "Dim1", y = "Dim2")
ggsave(file.path(out_dir, "R_embedding.png"), p, width = 6, height = 4, dpi = 100)

# Imbalance score
scores <- imbalance_score(Object = as.matrix(toy$sd[, c("Dim1", "Dim2")]),
                          conditions = as.character(toy$sd$conditions))
df$score <- scores$scaled_scores
p2 <- ggplot(df, aes(x = Dim1, y = Dim2, col = score)) +
  geom_point(size = 1.5) +
  theme_classic() +
  scale_color_viridis_c()
ggsave(file.path(out_dir, "R_imbalance.png"), p2, width = 6, height = 4, dpi = 100)

# Save data for Py
write.csv(df, file.path(out_dir, "df.csv"), row.names = FALSE)
write.csv(toy$sd, file.path(out_dir, "sd.csv"), row.names = FALSE)
cat("R plots done\n")
