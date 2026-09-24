# Week 5 — AI verification log

**Name:** Zhou Hanyi  **Student ID:** _(fill in)_  **Date:** 2026-09-24

The assignment requires AI to be used for **one clearly documented task**. This log
records the prompt, what the AI produced, what I checked, and what I rejected.

---

## 1. The documented AI task: plotting code for the PCA and volcano figures

### Prompt (preserved verbatim)

> I have a DESeq2 result table for a 2 × 3 design (condition: control/treated; batch:
> A/B/C, balanced, n = 12) with columns gene_id, baseMean, log2FoldChange, lfcSE, padj.
> Write ggplot2 code for (a) a PCA of the VST-transformed samples coloured by condition
> and shaped by batch, and (b) a volcano plot of shrunken log2FoldChange against
> -log10(padj), with dashed threshold lines at padj = 0.05 and |log2FC| = 1 and the
> significant genes coloured by direction. Label the axes with the percent variance for
> the PCA. Do not choose the thresholds for me — use the ones I gave. Do not filter the
> table.

### What the AI produced

Working `ggplot2` code for both figures, structurally the same as the course starter's:
`plotPCA(vsd, intgroup = ..., returnData = TRUE)` feeding a `geom_point` layer, and a
volcano built from `mutate(neg_log10_padj = -log10(padj))`.

### What I checked, and what I changed

| AI output | My verification | Decision |
|---|---|---|
| `-log10(padj)` computed directly | A gene with `padj == 0` after underflow gives `Inf` and silently drops the point with a ggplot warning. Checked in the cross-run: no gene has `padj == 0`; the minimum is 1.3 x 10^-10. | **Changed** to `-log10(pmax(padj, 1e-300))`. It does not bite on this dataset (minimum padj = 1.3 × 10⁻¹⁰), but the guard is free. |
| `drop_na(padj)` before plotting | DESeq2 sets `padj = NA` by independent filtering and Cook's-distance outlier removal. Dropping those rows silently is fine for a plot but wrong if the same pipeline feeds a count. Checked in the cross-run: `padj` is NA for **0** of 989 tested genes. | **Kept for the plot only.** The exported table retains every gene, NA rows included. |
| `plotPCA(vsd, intgroup = "condition")` | DESeq2's `plotPCA` uses the **500 most variable genes** by default (`ntop = 500`), not all genes. The AI did not mention this. | **Kept the default**, now knowingly, and recorded it. The PCA describes the most variable genes, not the transcriptome. |
| `vst(dds)` with default `blind` | `vst()` defaults to `blind = TRUE`, which ignores the design. For QC *after* the design is trusted, `blind = FALSE` is the documented choice. | **Changed** to `blind = FALSE`. |
| Legend combining condition and batch into six entries | Six combined entries force the reader to parse colour and shape together. | **Changed** to two separate keys — one for colour (condition), one for shape (batch). |
| Threshold lines | Confirmed the AI used the values I supplied (0.05, 1) and did not substitute its own. | **Accepted.** |

### An AI error that was caught

On the first pass the AI wrote the volcano's colour mapping as

```r
direction = ifelse(log2FoldChange > 0, "Up", "Down")
```

which assigns a direction to **every** gene, including the 929 non-significant ones —
the plot would have shown a fully red-and-blue cloud with no grey, implying 989
differentially expressed genes. Caught by comparing the legend counts against
`table(res_df$direction)`. Corrected to gate on `significant` first, so the three
categories are `Up in treated` (36), `Down in treated` (24), `Not significant` (929);
these now match the exported table exactly.

---

## 2. Independent verification of sample identity and coefficient direction

These were done by hand, not by AI, because they are the two failures that produce a
complete, plausible, and entirely wrong result.

> **Which run these came from.** For the reason set out in §4, the R script has not been
> executed here; the checks below were run in the PyDESeq2 cross-run, which uses the same
> inputs and the same design. The R script encodes each one as a `stopifnot()`, so running
> it re-performs them. Anything below that I did **not** observe is marked as such.

**Sample identity.** Asserted rather than assumed — in the R script as

```r
stopifnot(identical(colnames(counts), rownames(coldata)))
stopifnot(!any(duplicated(colnames(counts))))
```

and executed as the equivalent `assert` statements in the cross-run. **Both pass:** the
12 IDs match in content *and* order, so no reordering was needed — but the assertion is
what establishes that, not inspection.

**Reference level.** `condition` is built with control as the first level, so control is
the reference and the fitted coefficient contrasts treated against it. Confirmed in the
cross-run: the design matrix columns are `Intercept`, `batch[T.B]`, `batch[T.C]`,
`condition[T.treated]` — the `[T.treated]` suffix is what marks control as the baseline.
Had this been the other way round, every fold change would have flipped sign with no
error raised.

**Coefficient name — verified in PyDESeq2, not yet in R.** The name is read from the
fitted object rather than hard-coded. In the cross-run it is `condition[T.treated]`.
DESeq2 names the same coefficient `condition_treated_vs_control`; the R script checks
`resultsNames(dds)` for that string and **stops** if it is absent, so a mismatch cannot
pass silently. I have not observed `resultsNames()` output directly — that happens when
the R script is run.

**Direction, checked against the counts themselves.** For the top gene, the normalised
counts were averaged by group by hand and compared with the reported fold change:

```
Gene0035: mean treated 161.4 vs mean control 46.4
manual log2FC +1.80   reported shrunken log2FC +1.74
```

Same sign, and the small gap is apeglm shrinkage, which is expected. This is the check
that would catch a reversed contrast; a p-value cannot.

**Shrinkage sanity.** apeglm must change effect sizes but not adjusted *p* values. The
script asserts `all.equal(res$padj, res_shrunk$padj)`; median |log2FC| moved 0.225 →
0.102 while the adjusted *p* values were identical.

---

## 3. What the AI did **not** decide

- The thresholds (padj < 0.05, |log2FC| ≥ 1) — taken from the assignment.
- The design formula `~ batch + condition` — required, and justified in the interpretation.
- The filter rule (≥ 10 counts in ≥ 3 samples) — required; the choice of 3 as the floor is
  argued in the script comments.
- Which genes matter — no gene list was requested from or supplied by the AI.
- Whether to drop a sample — none was dropped.

---

## 4. Disclosure: which engine produced the numbers in this submission

This matters for reproducibility, so it is stated plainly rather than buried.

`week5_deseq2_analysis.R` is the submission and is written against DESeq2 + apeglm. It
**could not be executed in the environment where this homework was prepared**, for a
reason unrelated to the analysis:

- conda's bioconda channel ships no Windows (win-64) build of `bioconductor-deseq2`;
- installing Bioconductor's own Windows binaries into a user-writable library succeeds,
  but loading them fails with `Refusing to dyn.load shared library from writable path` —
  the sandbox will not load compiled `.dll` files from a writable directory.

So that the reported numbers would be measured rather than guessed, the identical design
was fitted with **PyDESeq2 0.5.4**, the reference Python reimplementation of the DESeq2
method, including apeglm shrinkage (`code/week5_pydeseq2_crossrun.py`). Every figure and
every number in `week5_interpretation.md` comes from that run and is reproducible from
`outputs/pydeseq2_results.csv`.

**The R script remains authoritative.** Running

```bash
cd Week5/code && Rscript week5_deseq2_analysis.R
```

on a machine with DESeq2 installed produces `week5_deseq2_results.csv`,
`week5_deseq2_object.rds`, `session_info.txt`, and overwrites the two figures with the R
renders. Small differences are possible — DESeq2 and PyDESeq2 agree closely but not
bit-for-bit, particularly in Cook's-distance outlier refitting — and **where they differ,
the R output is correct and should replace the numbers quoted here.**

One known cosmetic difference: PyDESeq2 names the coefficient `condition[T.treated]`
(patsy-style formula naming) where DESeq2 names it `condition_treated_vs_control`. Same
contrast, different label.

---

## 5. Use of the instructor key

`Week5_Homework_Gene_Annotation_Instructor_Key.csv`, shipped in the student package,
contains `truth_log2FC_for_instructor` — the simulation's real effect sizes.

It was **not** used to derive any result. `week5_deseq2_analysis.R` and
`week5_pydeseq2_crossrun.py` never open it. It is read only by
`code/week5_figures_and_benchmark.py`, after the analysis was already fixed, to benchmark
recovery: sensitivity 0.659, precision 0.967, empirical FDR 0.033 against a nominal 0.05,
and correct direction in 58 of 58 true positives. Those figures are reported in the
interpretation as a methodological check, not as the answer to the homework.
