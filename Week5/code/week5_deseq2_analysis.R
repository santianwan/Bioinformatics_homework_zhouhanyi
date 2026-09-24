# =============================================================================
# Week 5 Homework 1 — From verified counts to an interpretable DESeq2 result
# Bioinformatics: From Multi-Omics Data to Discovery
#
# Built from Week5_Homework_Starter.R. Differences from the starter are marked
# [CHANGED] with the reason; everything else follows the required workflow.
#
# Run from Week5/code/ :   Rscript week5_deseq2_analysis.R
# Inputs  : ../data/Week5_Homework_Count_Matrix.csv
#           ../data/Week5_Homework_Sample_Metadata.csv
# Outputs : ../outputs/  and  ../figures/
#
# This script is BLIND to Week5_Homework_Gene_Annotation_Instructor_Key.csv.
# That file carries the simulation's true effect sizes; it is read only by the
# separate script week5_truth_benchmark.R, AFTER this analysis is fixed.
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(tidyverse)
  library(ggrepel)
})

# [CHANGED] The starter uses rstudioapi::getActiveDocumentContext() to set the
# working directory. That only works inside an interactive RStudio session and
# errors under Rscript, so paths are relative to this script's location instead.
count_file    <- file.path("..", "data", "Week5_Homework_Count_Matrix.csv")
metadata_file <- file.path("..", "data", "Week5_Homework_Sample_Metadata.csv")
out_dir <- file.path("..", "outputs")
fig_dir <- file.path("..", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(1)   # nothing here is stochastic, but fixed for reproducibility

# ----------------------------------------------------------------- 1. Import
counts_df <- read.csv(count_file, row.names = 1, check.names = FALSE)
coldata   <- read.csv(metadata_file, row.names = 1, check.names = FALSE)
counts    <- as.matrix(counts_df)

# ------------------------------------------------------- 2. Input validation
# Each check is a separate stopifnot() so a failure names the assumption that
# broke, rather than reporting "one of four things is wrong".
stopifnot("column count != sample count" = ncol(counts) == nrow(coldata))
stopifnot("sample IDs differ or are out of order" =
            identical(colnames(counts), rownames(coldata)))
stopifnot("duplicated sample IDs" = !any(duplicated(colnames(counts))))
stopifnot("duplicated gene IDs" = !any(duplicated(rownames(counts))))
stopifnot("missing values in the count matrix" = !any(is.na(counts)))
stopifnot("negative counts" = all(counts >= 0))
stopifnot("non-integer counts" = all(counts == round(counts)))

# [CHANGED] Added an explicit guard that the matrix is raw counts, not a
# normalised or transformed matrix. A TPM/CPM/z-scored matrix can still pass
# "non-negative" but would invalidate the negative-binomial model. Raw count
# columns do not sum to a constant and do contain zeros.
stopifnot("no zero counts anywhere - matrix may not be raw" = any(counts == 0))
lib <- colSums(counts)
stopifnot("column sums are suspiciously constant - matrix may be normalised" =
            (max(lib) - min(lib)) / mean(lib) > 1e-6)

coldata$condition <- relevel(factor(coldata$condition), ref = "control")
coldata$batch     <- factor(coldata$batch)
stopifnot("unexpected condition levels" =
            identical(levels(coldata$condition), c("control", "treated")))
stopifnot("unexpected batch levels" =
            identical(levels(coldata$batch), c("A", "B", "C")))

cat("=== Design ===\n");  print(table(coldata$batch, coldata$condition))
cat("\n=== Library sizes (total assigned counts per sample) ===\n")
print(summary(lib))
cat("range:", format(min(lib), big.mark = ","), "-",
    format(max(lib), big.mark = ","),
    sprintf("(%.2fx)\n", max(lib) / min(lib)))
cat("genes detected per sample:", min(colSums(counts > 0)), "-",
    max(colSums(counts > 0)), "of", nrow(counts), "\n\n")

# ------------------------------------------------ 3. Construct DESeq2 object
# WHY batch is in the design: the three batches are balanced (2 control +
# 2 treated each), so batch is not confounded with condition and the treatment
# estimate would remain unbiased without it. It is included because it absorbs
# a real source of between-sample variance: leaving it in the residual inflates
# the dispersion estimate and costs power. Balance protects the point estimate;
# modelling the batch protects the standard error.
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = coldata,
                              design    = ~ batch + condition)

# [CHANGED] Explicit full-rank check. With 3 batches + 2 conditions and 12
# samples the model matrix has 4 columns and must be rank 4; a rank-deficient
# matrix is the failure mode that silently follows an unbalanced design.
mm <- model.matrix(design(dds), data = as.data.frame(colData(dds)))
stopifnot("model matrix is rank deficient" = qr(mm)$rank == ncol(mm))
cat("model matrix:", nrow(mm), "samples x", ncol(mm), "coefficients, rank",
    qr(mm)$rank, "(full rank)\n\n")

# --------------------------------------------------------------- 4. Filtering
# Rule: keep genes with at least 10 counts in at least 3 samples.
# 3 is the floor that a gene expressed in only one batch-by-condition cell
# (n = 2) cannot reach, so single-cell-of-the-design artefacts are removed
# while a gene present in any one full condition group (n = 6) is retained.
keep <- rowSums(counts(dds) >= 10) >= 3
n_before <- nrow(dds)
dds <- dds[keep, ]
cat("=== Filtering: >= 10 counts in >= 3 samples ===\n")
cat("genes before:", n_before, "\ngenes after: ", nrow(dds),
    sprintf("(%.1f%% removed)\n\n", 100 * (1 - nrow(dds) / n_before)))

# -------------------------------------------------------------- 5. Fit model
dds <- DESeq(dds)
coef_names <- resultsNames(dds)
cat("=== resultsNames(dds) ===\n"); print(coef_names); cat("\n")

# The coefficient name is READ from the fitted object, not assumed. It is only
# checked against the expected string so that a silent change in level ordering
# (which would flip the sign of every fold change) stops the script.
target_coef <- "condition_treated_vs_control"
if (!target_coef %in% coef_names) {
  stop("Expected coefficient not found. resultsNames(dds) = ",
       paste(coef_names, collapse = ", "))
}

# ------------------------------------------- 6/7. Extract contrast, then shrink
res <- results(dds, contrast = c("condition", "treated", "control"), alpha = 0.05)
cat("=== Unshrunken result ===\n"); print(summary(res))

res_shrunk <- lfcShrink(dds, coef = target_coef, type = "apeglm")

# [CHANGED] Verify that shrinkage changed only the effect sizes, not the
# inference. apeglm keeps the adjusted p values of results(); if these differ,
# the wrong coefficient or contrast was shrunk.
stopifnot("apeglm changed the adjusted p values - wrong coefficient?" =
            isTRUE(all.equal(res$padj, res_shrunk$padj)))
cat("\nshrinkage: median |log2FC| ",
    sprintf("%.3f -> %.3f", median(abs(res$log2FoldChange), na.rm = TRUE),
            median(abs(res_shrunk$log2FoldChange), na.rm = TRUE)), "\n")

# [CHANGED] Direction sanity check, independent of DESeq2's labelling. For the
# top gene, compare the mean normalised count in treated vs control by hand and
# confirm the sign agrees with the reported log2FoldChange. This is the check
# that catches a reversed contrast, which no amount of p-value inspection will.
norm_counts <- counts(dds, normalized = TRUE)
top_gene <- rownames(res_shrunk)[which.min(res_shrunk$padj)]
mean_t <- mean(norm_counts[top_gene, coldata$condition == "treated"])
mean_c <- mean(norm_counts[top_gene, coldata$condition == "control"])
manual_lfc <- log2(mean_t / mean_c)
cat(sprintf("direction check on %s: mean treated %.1f vs control %.1f, ",
            top_gene, mean_t, mean_c),
    sprintf("manual log2FC %+.2f, DESeq2 shrunken %+.2f\n",
            manual_lfc, res_shrunk[top_gene, "log2FoldChange"]))
stopifnot("reported fold-change sign disagrees with the normalised counts" =
            sign(manual_lfc) == sign(res_shrunk[top_gene, "log2FoldChange"]))

# --------------------------------------------------------------- 8. Result table
PADJ_CUT <- 0.05
LFC_CUT  <- 1

res_df <- as.data.frame(res_shrunk) |>
  rownames_to_column("gene_id") |>
  mutate(
    significant = !is.na(padj) & padj < PADJ_CUT & abs(log2FoldChange) >= LFC_CUT,
    direction = case_when(
      significant & log2FoldChange > 0 ~ "Up in treated",
      significant & log2FoldChange < 0 ~ "Down in treated",
      TRUE                             ~ "Not significant"
    )
  ) |>
  arrange(padj)

# The FULL table is exported, non-significant genes included: a results file
# that contains only the hits cannot be used to check anything.
write.csv(res_df, file.path(out_dir, "week5_deseq2_results.csv"), row.names = FALSE)

cat("\n=== Significance: padj <", PADJ_CUT, "AND |log2FC| >=", LFC_CUT, "===\n")
print(table(res_df$direction))
cat("padj < 0.05 regardless of effect size:",
    sum(!is.na(res_df$padj) & res_df$padj < PADJ_CUT), "\n")
cat("padj = NA (independent filtering):", sum(is.na(res_df$padj)), "\n\n")

# ----------------------------------------------------------------- 9. PCA
# blind = FALSE: the design is already known and trusted here, so the
# variance-stabilising fit may use it. Transformed values are for visualisation
# only and are never fed back into testing.
vsd <- vst(dds, blind = FALSE)
pca_df <- plotPCA(vsd, intgroup = c("condition", "batch"), returnData = TRUE)
percent_var <- round(100 * attr(pca_df, "percentVar"))

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition, shape = batch, label = name)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3, max.overlaps = Inf, show.legend = FALSE) +
  scale_color_manual(values = c(control = "#2F6DB3", treated = "#C0392B")) +
  labs(title = "Condition separates on PC1; batch does not structure the samples",
       subtitle = "Variance-stabilising transform, 12 samples, balanced 2x3 design",
       x = paste0("PC1: ", percent_var[1], "% variance"),
       y = paste0("PC2: ", percent_var[2], "% variance"),
       color = "Condition", shape = "Batch") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(size = 11), legend.position = "right")

ggsave(file.path(fig_dir, "week5_pca.png"), p_pca, width = 7.2, height = 5, dpi = 300)

# Quantify what the PCA shows, rather than asserting it from the picture.
cat("=== PC1 by group ===\n")
print(pca_df |> group_by(condition) |>
        summarise(mean_PC1 = mean(PC1), .groups = "drop"))
print(pca_df |> group_by(batch) |>
        summarise(mean_PC1 = mean(PC1), mean_PC2 = mean(PC2), .groups = "drop"))
cat("\nPC1 ~ condition:", format.pval(
  summary(aov(PC1 ~ condition, data = pca_df))[[1]][["Pr(>F)"]][1], digits = 3), "\n")
cat("PC1 ~ batch:    ", format.pval(
  summary(aov(PC1 ~ batch, data = pca_df))[[1]][["Pr(>F)"]][1], digits = 3), "\n")
cat("PC2 ~ batch:    ", format.pval(
  summary(aov(PC2 ~ batch, data = pca_df))[[1]][["Pr(>F)"]][1], digits = 3), "\n\n")

# --------------------------------------------------------------- 10. Volcano
plot_df <- res_df |> mutate(neg_log10_padj = -log10(pmax(padj, 1e-300)))
label_df <- plot_df |> filter(significant) |> slice_head(n = 10)

p_volcano <- ggplot(plot_df, aes(log2FoldChange, neg_log10_padj, color = direction)) +
  geom_point(alpha = 0.7, size = 1.6) +
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = -log10(PADJ_CUT), linetype = "dashed", colour = "grey40") +
  geom_text_repel(data = label_df, aes(label = gene_id), size = 2.6,
                  max.overlaps = Inf, show.legend = FALSE, colour = "black") +
  annotate("text", x = Inf, y = -log10(PADJ_CUT), hjust = 1.05, vjust = -0.6,
           label = "padj = 0.05", size = 2.8, colour = "grey30") +
  annotate("text", x = LFC_CUT, y = Inf, hjust = -0.15, vjust = 1.6,
           label = "|log2FC| = 1", size = 2.8, colour = "grey30") +
  scale_color_manual(values = c("Up in treated"    = "#C0392B",
                                "Down in treated"  = "#2F6DB3",
                                "Not significant"  = "grey75")) +
  labs(title = "Treated versus control, after apeglm shrinkage",
       subtitle = paste0("Thresholds: adjusted p < ", PADJ_CUT,
                         " and |log2 fold change| >= ", LFC_CUT,
                         ";  ", sum(res_df$significant), " of ", nrow(res_df),
                         " tested genes significant"),
       x = "Shrunken log2 fold change (treated / control)",
       y = "-log10 adjusted p value", color = NULL) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(size = 11), legend.position = "bottom")

ggsave(file.path(fig_dir, "week5_de_plot.png"), p_volcano, width = 7.2, height = 5.4, dpi = 300)

# ------------------------------------------------------- 11. Reproducibility
saveRDS(dds, file.path(out_dir, "week5_deseq2_object.rds"))
capture.output(sessionInfo(), file = file.path(out_dir, "session_info.txt"))

# A compact machine-readable summary so the interpretation can quote numbers
# without anyone re-deriving them by hand.
summary_list <- list(
  genes_input = n_before, genes_tested = nrow(dds),
  samples = ncol(dds), design = "~ batch + condition",
  coefficient = target_coef,
  padj_cut = PADJ_CUT, lfc_cut = LFC_CUT,
  n_significant = sum(res_df$significant),
  n_up = sum(res_df$direction == "Up in treated"),
  n_down = sum(res_df$direction == "Down in treated"),
  n_padj_only = sum(!is.na(res_df$padj) & res_df$padj < PADJ_CUT),
  n_padj_na = sum(is.na(res_df$padj)),
  lib_min = min(lib), lib_max = max(lib),
  pc1_var = percent_var[1], pc2_var = percent_var[2],
  top_gene = top_gene
)
writeLines(jsonlite::toJSON(summary_list, auto_unbox = TRUE, pretty = TRUE),
           file.path(out_dir, "week5_run_summary.json"))

cat("=== Wrote ===\n")
cat(" outputs/week5_deseq2_results.csv  (", nrow(res_df), " genes, full table )\n", sep = "")
cat(" outputs/week5_deseq2_object.rds\n outputs/session_info.txt\n")
cat(" outputs/week5_run_summary.json\n")
cat(" figures/week5_pca.png\n figures/week5_de_plot.png\n")
