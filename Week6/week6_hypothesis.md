# Week 6 — 16S Microbiome Analysis in EasyMultiProfiler-Web

**Name:** Zhou Hanyi  **Student ID:** SUAT24000202  **Date:** 2026-09-25

Homework (from `Homework-Week6.docx`): use the 16S files in the EMP `tests` folder,
run the full analysis in EasyMultiProfiler-Web, submit via Sync, and propose a
scientific hypothesis from the results, naming which parameters back it up.

## What's in this dataset

`16S_level-7.csv` + `16S_mapping.csv`: 130 stool samples, 470 genus/species-level
taxa. The `Group` column packs three things into one label —
`{IBS,UC}_{before,after}` — and `Group_sub` adds a `great`/`poor` treatment-response
tag on top of that. So really it's a 2×2×2 design: disease × timepoint × response,
roughly balanced (13–18 samples per cell — UC has fewer patients than IBS, which
matters later).

I didn't want to walk into the web tool guessing what the numbers should look like,
so I fit the same table independently first with `vegan` in R — Shannon diversity,
Bray-Curtis + PERMANOVA, and a Wilcoxon screen for baseline taxa. Script and full
output below. That's what the "which parameters support your hypothesis" part of the
report is actually built on.

## What the independent run found

Kruskal-Wallis on Shannon diversity: disease matters (p = 0.019), timepoint alone
doesn't (p = 0.22), response alone doesn't (p = 0.10). PERMANOVA on Bray-Curtis
tells the same story from the composition side — disease is the only term that
clears p < 0.05 (R² = 2.9%, p = 0.001), and every interaction term, including
disease×timepoint×response, comes out non-significant.

That's a smaller effect than I expected going in — I'd guessed timepoint (treatment)
would be the dominant axis, since that's usually the headline in a before/after
design. It isn't, at least not overall. Splitting by disease instead:

- **Within UC**, Shannon drops from before to after (means 2.77 → 2.43 across
  response groups), and the test is close to conventional significance
  (p = 0.057, n = 58). Not clean enough to call confirmed on its own.
- **Within IBS**, there's no such drop at all (p = 0.80, n = 72) — before and after
  are essentially the same distribution.

So whatever timepoint is doing, it isn't doing it evenly across the two diseases,
and UC is the one where it looks real. The boxplot below shows this directly: the
UC-after box sits visibly lower than the other three, including outliers down near
Shannon = 1.

None of the 155 taxa I tested for a before-treatment poor-vs-great difference
survived BH correction (best padj was 1.0 across the board, despite raw p-values as
low as 0.02). A few names came up worth keeping in mind if a larger cohort ever
tests this again — *Akkermansia muciniphila* and *Lactobacillus salivarius* trended
higher in great responders, both are commonly reported as beneficial commensals —
but at n≈65 per timepoint this is nowhere near enough power to call it.

![PCoA](figures/week6_pcoa.png)

*Bray-Curtis PCoA, coloured by disease, shaped by timepoint. IBS and UC overlap
substantially but UC pulls toward higher PCoA1 — matches the PERMANOVA result of a
real but modest (R²=2.9%) disease effect, not a clean two-cluster split.*

![Alpha diversity](figures/week6_alpha_boxplot.png)

*Shannon diversity by disease and timepoint. UC starts lower than IBS and drops
further after treatment; IBS barely moves.*

## Checking that against what EasyMultiProfiler-Web actually produced

I didn't want to just assume the web tool would reproduce this, so after running its
One-click pipeline (genus level, Shannon, Bray-Curtis/PCoA — same choices as above) I
pulled its own `02_alpha_indices.csv` back out of the exported bundle and reran the
identical Kruskal-Wallis tests on the web tool's own Shannon numbers, not mine. Script:
[`code/week6_webtool_crosscheck.R`](code/week6_webtool_crosscheck.R).

This did not come out the way I expected. Per-sample, the two Shannon calculations
agree closely in rank (Spearman ρ = 0.81 across 130 matched samples, p < 2.2×10⁻¹⁶) —
so it's the same underlying signal, not a different computation. But the group
comparisons that were significant in my untouched run are not significant on the web
tool's own numbers:

| Test | My reference (raw counts, no filter) | Web tool's own numbers |
|---|---:|---:|
| Shannon ~ disease | p = 0.019 | p = 0.152 |
| Shannon ~ timepoint, within UC | p = 0.057 | p = 0.211 |

The web tool's Preprocess step applied its own default low-abundance filter (min
prevalence 0.1, min detect rate 0.05) plus a genus-level taxonomy collapse before
computing diversity — steps I deliberately skipped in the independent check so I'd
have an unfiltered baseline to compare against. That combination visibly changes
per-sample richness (e.g. sample J_XYL_F_0001_01: 25 raw observed genera/species drop
to 8 after the tool's own pipeline) without reversing the direction of any group
difference — the UC-after box is still visibly the lowest of the four in the tool's
own boxplot (median ≈1.44 vs 1.58–1.65 for the other three; see
`outputs/webtool_run_20260925-024344/plots/02_alpha_shannon_boxplot.png`) — but it
does drop two of my three key comparisons below conventional significance.

So there are really two findings here, not one, and the second one wasn't the one I
went in looking for: the disease/diversity relationship itself, and the fact that a
default preprocessing filter is not a neutral cleaning step — it can erase a real
between-group signal by removing exactly the low-abundance taxa that were carrying
it. That is worth reporting as an actual result, not folded quietly into a methods
paragraph.

## The hypothesis

**Ulcerative colitis carries a baseline gut-microbiome diversity deficit relative to
IBS, and this deficit does not recover with treatment — if anything it may deepen —
but detecting it is sensitive to preprocessing choices that are easy to apply without
thinking about them.** On the raw genus/species table, Shannon diversity differs by
disease (p = 0.019), and within UC alone keeps falling from before to after treatment
(p = 0.057) while IBS shows no such trend (p = 0.80); PERMANOVA on the same raw table
finds disease as the one significant driver of Bray-Curtis dissimilarity (R² = 2.9%,
p = 0.001). Running the identical comparisons on EasyMultiProfiler-Web's own output —
after its default prevalence/detect-rate filter and genus-level collapse — keeps the
same direction in every case but pushes both p-values above 0.05. The taxa doing the
filtering-out are apparently informative, not just noise.

I'd stop short of calling the UC before→after drop confirmed on either version of the
analysis — 0.057 was never below 0.05 to begin with, and UC has the smaller n of the
two diseases (58 vs 72 samples), which cuts both ways: real effects are harder to
detect there, but so is separating a real trend from noise. What I'm more confident
about is the disease effect itself, since it shows up two different ways (alpha and
beta diversity) on the raw data, and the absence of any detectable baseline
response-group signal, on either version of the pipeline — good/poor responders look
the same on 16S alone before treatment starts.

If I were extending this: rerunning the web tool's pipeline with a looser prevalence
threshold (or none at all) would directly test whether the filter is really what's
costing the significance, rather than inferring it the way I did here from a rank
correlation.

---

## What was done in EasyMultiProfiler-Web

- **Import:** "16S Microbiome (Course Demo)" one-click loader — confirmed against the
  tool's own `demo_data.R` that this reads `tests/16S_level-7.csv` and
  `tests/16S_mapping.csv` directly, the same files named in the assignment. Loaded as
  132 samples × 470 features (2 of the 132 assay columns don't have a mapping-file
  match — a data quirk, not something I introduced).
- **Preprocess:** filter tab, default values — MIN MAX COUNT 0, MIN DETECT RATE 0.05,
  MAX DETECT RATE 1, MIN PREVALENCE 0.1, MAX NA PROPORTION 1.
- **Analysis:** the microbiome One-click pipeline — genus-level, Shannon alpha
  diversity, Bray-Curtis beta diversity with PCoA/PCA/NMDS ordination, |log2FC| ≥ 1
  and p ≤ 0.05 for the (empty-output) differential-taxa step.
- **Sync:** submitted to GitHub — commit and run path below.

Full output from this run (summary log, alpha/beta tables, plots) is kept in
[`outputs/webtool_run_20260925-024344/`](outputs/webtool_run_20260925-024344/).

## Submission record

- **Synced to GitHub:** commit
  [`4e98dea`](https://github.com/santianwan/Bioinformatics_homework_zhouhanyi/commit/4e98dea72b70379cc139779795d0a06be0c58092),
  landed at `EMP2026/Week_06/microbiome_16s/weekly/runs/2026-09-24T18-53-27-018Z-11jyf8/`.
  Confirmed the synced `results/m16s_course_alpha.csv` matches the locally downloaded
  ZIP bundle exactly (same per-sample Shannon values, e.g. J_XYL_F_0001_01 = 0.9231 in
  both) — it's the same run, not a different one overwriting it.
- The first Sync attempt failed with `Could not snapshot session file
  empt_rnaseq_course.rds ... embedded nul in string` — a leftover, corrupted session
  snapshot from the unrelated Week 5 RNA-seq work still sitting in the same EMP-Web
  session. Clearing all loaded experiments and reloading only the 16S Course Demo
  before re-running Preprocess and One-click Run resolved it. (Week 5's RNA-seq result
  was already synced separately before this, so clearing it here lost nothing.)

---

## Appendix: independent reference analysis script

```r
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

write.csv(res_tab, "outputs/week6_baseline_response_taxa.csv", row.names = FALSE)
write.csv(meta[, c("Group","Group_sub","disease","timepoint","response","shannon","observed")],
          "outputs/week6_alpha_diversity.csv")

# ---- figures ----
library(ggplot2)

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
ggsave("figures/week6_pcoa.png", p_pcoa, width = 6.5, height = 4.6, dpi = 300)

p_alpha <- ggplot(meta, aes(disease, shannon, fill = timepoint)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.85) +
  labs(title = "UC carries lower Shannon diversity than IBS at both timepoints",
       subtitle = "Kruskal-Wallis: disease p = 0.019; timepoint p = 0.22 (ns overall)",
       x = NULL, y = "Shannon diversity", fill = "Timepoint") +
  theme_bw(base_size = 11)
ggsave("figures/week6_alpha_boxplot.png", p_alpha, width = 6, height = 4.4, dpi = 300)
```
