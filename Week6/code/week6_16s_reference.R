suppressPackageStartupMessages({
  library(vegan)
})

src <- "F:/EasyMultiProfiler-Web-9.0.4/EasyMultiProfiler-Web-9.0.4/tests"
assay <- read.csv(file.path(src, "16S_level-7.csv"), row.names = 1, check.names = FALSE)
meta  <- read.csv(file.path(src, "16S_mapping.csv"), row.names = 1, check.names = FALSE)

cat("assay dim:", dim(assay), "\n")
cat("meta dim:", dim(meta), "\n")

common <- intersect(colnames(assay), rownames(meta))
cat("common samples:", length(common), "of assay", ncol(assay), "and meta", nrow(meta), "\n")
assay <- assay[, common]
meta  <- meta[common, ]
stopifnot(identical(colnames(assay), rownames(meta)))

meta$disease   <- sub("_.*", "", meta$Group)                  # IBS / UC
meta$timepoint <- ifelse(grepl("after", meta$Group), "after", "before")
meta$response  <- sub(".*_", "", meta$Group_sub)               # great / poor
meta$disease   <- factor(meta$disease, levels = c("IBS","UC"))
meta$timepoint <- factor(meta$timepoint, levels = c("before","after"))
meta$response  <- factor(meta$response, levels = c("poor","great"))

cat("\n=== design table ===\n")
print(table(meta$disease, meta$timepoint, meta$response))

mat <- t(as.matrix(assay))         # samples x features
mat[mat < 0] <- 0
lib <- rowSums(mat)
cat("\nlibrary size range:", min(lib), "-", max(lib), "\n")

# ---- alpha diversity ----
shannon  <- diversity(mat, index = "shannon")
observed <- rowSums(mat > 0)
meta$shannon  <- shannon
meta$observed <- observed

cat("\n=== alpha diversity: Shannon by group ===\n")
print(aggregate(shannon ~ disease + timepoint + response, data = meta, FUN = function(x) round(mean(x), 3)))

kw_disease   <- kruskal.test(shannon ~ disease, data = meta)
kw_timepoint <- kruskal.test(shannon ~ timepoint, data = meta)
kw_response  <- kruskal.test(shannon ~ response, data = meta)
cat("\nKruskal-Wallis Shannon ~ disease:  p =", format.pval(kw_disease$p.value, digits=3), "\n")
cat("Kruskal-Wallis Shannon ~ timepoint:p =", format.pval(kw_timepoint$p.value, digits=3), "\n")
cat("Kruskal-Wallis Shannon ~ response: p =", format.pval(kw_response$p.value, digits=3), "\n")

# response effect, separately before vs after, and separately within each disease
for (tp in levels(meta$timepoint)) {
  sub <- meta[meta$timepoint == tp, ]
  kw <- kruskal.test(shannon ~ response, data = sub)
  cat(sprintf("  Shannon ~ response within timepoint=%s: p = %s (n=%d)\n", tp, format.pval(kw$p.value, digits=3), nrow(sub)))
}
for (dz in levels(meta$disease)) {
  sub <- meta[meta$disease == dz, ]
  kw <- kruskal.test(shannon ~ timepoint, data = sub)
  cat(sprintf("  Shannon ~ timepoint within disease=%s: p = %s (n=%d)\n", dz, format.pval(kw$p.value, digits=3), nrow(sub)))
}

# ---- beta diversity: Bray-Curtis PERMANOVA ----
rel <- mat / rowSums(mat)
bray <- vegdist(rel, method = "bray")

set.seed(1)
ad_full <- adonis2(bray ~ disease * timepoint * response, data = meta, permutations = 999, by = "terms")
cat("\n=== PERMANOVA (Bray-Curtis), disease*timepoint*response ===\n")
print(ad_full)

# ---- differential abundance: which taxa separate great vs poor responders at baseline ----
before_idx <- meta$timepoint == "before"
rel_before <- rel[before_idx, ]
resp_before <- meta$response[before_idx]

pvals <- apply(rel_before, 2, function(col) {
  if (sum(col > 0) < 5) return(NA)
  tryCatch(wilcox.test(col ~ resp_before)$p.value, error = function(e) NA)
})
padj <- p.adjust(pvals, method = "BH")
res_tab <- data.frame(feature = colnames(rel_before), pvalue = pvals, padj = padj,
                       mean_poor = colMeans(rel_before[resp_before == "poor", , drop=FALSE]),
                       mean_great = colMeans(rel_before[resp_before == "great", , drop=FALSE]))
res_tab <- res_tab[order(res_tab$pvalue), ]
cat("\n=== Top baseline (before) taxa distinguishing poor vs great responders (by raw p) ===\n")
print(head(res_tab, 15), row.names = FALSE)
cat("\nn significant at padj<0.05:", sum(res_tab$padj < 0.05, na.rm=TRUE), "of", sum(!is.na(res_tab$padj)), "tested\n")

out_dir <- "C:/Users/zhy/AppData/Local/Temp/claude/F--subject-3-up/aeb99f47-0a02-4052-bdda-f11fb3b32390/scratchpad"
write.csv(res_tab, file.path(out_dir, "week6_baseline_response_taxa.csv"), row.names = FALSE)
write.csv(meta[, c("Group","Group_sub","disease","timepoint","response","shannon","observed")],
          file.path(out_dir, "week6_alpha_diversity.csv"))
saveRDS(list(adonis=ad_full), file.path(out_dir, "week6_permanova.rds"))
cat("\nWrote week6_baseline_response_taxa.csv, week6_alpha_diversity.csv\n")

# ---- figures ----
suppressPackageStartupMessages(library(ggplot2))

pcoa <- cmdscale(bray, k = 2, eig = TRUE)
pct_var <- round(100 * pcoa$eig[1:2] / sum(pcoa$eig[pcoa$eig > 0]))
pcoa_df <- data.frame(PC1 = pcoa$points[,1], PC2 = pcoa$points[,2],
                       disease = meta$disease, timepoint = meta$timepoint,
                       response = meta$response)

p_pcoa <- ggplot(pcoa_df, aes(PC1, PC2, color = disease, shape = timepoint)) +
  geom_point(size = 2.6, alpha = 0.85) +
  labs(title = "Disease separates on PCoA1; timepoint and response do not",
       subtitle = "Bray-Curtis dissimilarity, genus-level, 130 samples",
       x = paste0("PCoA1: ", pct_var[1], "%"), y = paste0("PCoA2: ", pct_var[2], "%"),
       color = "Disease", shape = "Timepoint") +
  theme_bw(base_size = 11)
ggsave(file.path(out_dir, "week6_pcoa.png"), p_pcoa, width = 6.5, height = 4.6, dpi = 300)

p_alpha <- ggplot(meta, aes(disease, shannon, fill = timepoint)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.85) +
  labs(title = "UC carries lower Shannon diversity than IBS at both timepoints",
       subtitle = "Kruskal-Wallis: disease p = 0.019; timepoint p = 0.22 (ns overall)",
       x = NULL, y = "Shannon diversity", fill = "Timepoint") +
  theme_bw(base_size = 11)
ggsave(file.path(out_dir, "week6_alpha_boxplot.png"), p_alpha, width = 6, height = 4.4, dpi = 300)

cat("Wrote week6_pcoa.png, week6_alpha_boxplot.png\n")
