# =============================================================================
# Week 5 — post-hoc benchmark against the simulation's true effect sizes
#
# IMPORTANT, and stated plainly: the student package ships
# `Week5_Homework_Gene_Annotation_Instructor_Key.csv`, whose column
# `truth_log2FC_for_instructor` is the ground-truth effect used to simulate the
# counts. That file is almost certainly meant for the instructor.
#
# It was NOT used to choose any threshold, filter, model term or gene list.
# week5_deseq2_analysis.R never opens it. This script runs strictly afterwards,
# on the already-written results table, and asks one question: how well did the
# workflow recover a signal whose true value is known?
#
# Run AFTER week5_deseq2_analysis.R, from Week5/code/ :
#   Rscript week5_truth_benchmark.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

res_path   <- file.path("..", "outputs", "week5_deseq2_results.csv")
truth_path <- file.path("..", "data", "Week5_Homework_Gene_Annotation_Instructor_Key.csv")
out_dir    <- file.path("..", "outputs")
fig_dir    <- file.path("..", "figures")

stopifnot("run week5_deseq2_analysis.R first" = file.exists(res_path))

res   <- read_csv(res_path, show_col_types = FALSE)
truth <- read_csv(truth_path, show_col_types = FALSE) |>
  rename(true_lfc = truth_log2FC_for_instructor)

# The key covers all 1,000 input genes; the results table covers only those that
# survived filtering. Join on the tested set and keep the filtered-out genes
# separately, because "missed because filtered" and "missed because not
# significant" are different failures.
bench <- res |> inner_join(truth |> select(gene_id, true_lfc), by = "gene_id")
filtered_out <- truth |> filter(!gene_id %in% res$gene_id)

stopifnot("join lost genes" = nrow(bench) == nrow(res))

# A gene counts as truly changed if its simulated effect clears the same
# effect-size threshold the analysis used. Using the analysis's own threshold
# keeps the comparison honest.
LFC_CUT <- 1
bench <- bench |>
  mutate(truly_changed = abs(true_lfc) >= LFC_CUT,
         called        = significant)

tp <- sum(bench$truly_changed & bench$called)
fp <- sum(!bench$truly_changed & bench$called)
fn <- sum(bench$truly_changed & !bench$called)
tn <- sum(!bench$truly_changed & !bench$called)

cat("=== Recovery of the simulated signal ===\n")
cat("tested genes:", nrow(bench), " (", nrow(filtered_out),
    "removed by the count filter before testing )\n", sep = "")
cat("truly changed (|true log2FC| >=", LFC_CUT, "):", sum(bench$truly_changed), "\n")
cat("called significant:", sum(bench$called), "\n\n")
cat(sprintf("  true positives  %4d\n  false positives %4d\n", tp, fp))
cat(sprintf("  false negatives %4d\n  true negatives  %4d\n\n", fn, tn))
cat(sprintf("sensitivity (recall)   %.3f\n", tp / (tp + fn)))
cat(sprintf("precision              %.3f\n", tp / (tp + fp)))
cat(sprintf("empirical FDR          %.3f   (nominal target 0.05)\n", fp / max(tp + fp, 1)))
cat(sprintf("specificity            %.3f\n\n", tn / (tn + fp)))

# Direction: of the genes correctly called, how many have the right sign?
dir_ok <- bench |> filter(called, truly_changed) |>
  summarise(n = n(), correct_sign = sum(sign(log2FoldChange) == sign(true_lfc)))
cat("sign of effect correct in", dir_ok$correct_sign, "of", dir_ok$n,
    "correctly called genes\n\n")

# Calibration of the shrunken estimate against the truth.
fit <- lm(log2FoldChange ~ true_lfc, data = bench)
cat("=== Calibration of shrunken log2FC against truth ===\n")
cat(sprintf("Pearson r      %.3f\n", cor(bench$log2FoldChange, bench$true_lfc)))
cat(sprintf("slope          %.3f  (1.0 = unbiased; <1 = shrinkage toward zero)\n",
            coef(fit)[2]))
cat(sprintf("intercept      %+.3f\n", coef(fit)[1]))
cat(sprintf("median |error| %.3f log2 units\n",
            median(abs(bench$log2FoldChange - bench$true_lfc))))
cat(sprintf("RMSE           %.3f\n\n", sqrt(mean((bench$log2FoldChange - bench$true_lfc)^2))))

# What did the count filter throw away?
if (nrow(filtered_out) > 0) {
  cat("=== Genes removed before testing ===\n")
  cat("removed:", nrow(filtered_out), "; of those, truly changed:",
      sum(abs(filtered_out$true_lfc) >= LFC_CUT), "\n")
  cat("median |true log2FC| among removed genes:",
      sprintf("%.3f\n\n", median(abs(filtered_out$true_lfc))))
}

# ------------------------------------------------------------------ figure
lab <- bench |> mutate(cls = case_when(
  truly_changed &  called ~ "Correctly called",
  truly_changed & !called ~ "Missed (false negative)",
 !truly_changed &  called ~ "False positive",
  TRUE                    ~ "Correctly not called"))

p <- ggplot(lab, aes(true_lfc, log2FoldChange, colour = cls)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45") +
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dotted", colour = "grey70") +
  geom_hline(yintercept = c(-LFC_CUT, LFC_CUT), linetype = "dotted", colour = "grey70") +
  geom_point(alpha = 0.75, size = 1.6) +
  scale_colour_manual(values = c("Correctly called"        = "#C0392B",
                                 "Missed (false negative)" = "#D98C1F",
                                 "False positive"          = "#6C3483",
                                 "Correctly not called"    = "grey75")) +
  labs(title = "Recovered effect sizes against the simulation's true values",
       subtitle = paste0("Dashed line is y = x (perfect recovery); fitted slope ",
                         sprintf("%.2f", coef(fit)[2]),
                         ", Pearson r ", sprintf("%.2f", cor(bench$log2FoldChange, bench$true_lfc))),
       x = "True log2 fold change (simulation)",
       y = "Shrunken log2 fold change (apeglm)", colour = NULL) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(size = 11), legend.position = "bottom")

ggsave(file.path(fig_dir, "week5_truth_benchmark.png"), p,
       width = 7.0, height = 5.4, dpi = 300)

write_csv(bench |> select(gene_id, baseMean, log2FoldChange, lfcSE, padj,
                          significant, direction, true_lfc, truly_changed),
          file.path(out_dir, "week5_truth_benchmark.csv"))

cat("=== Wrote ===\n outputs/week5_truth_benchmark.csv\n figures/week5_truth_benchmark.png\n")
