# Week 5 Homework 2 — Results and Interpretation (EasyMultiProfiler-Web)

**Name:** Zhou Hanyi  **Student ID:** SUAT24000202  **Date:** 2026-09-25

This is the write-up for [Week5_Homework_2_Instruction.md](https://github.com/xielab2017/Bioinformatics_SUAT_2026_FALL/blob/main/Week%205/Homework/for_student/Week5_Homework_2_Instruction.md):
run RNA-seq analysis in EasyMultiProfiler-Web, interpret the results, and Sync to GitHub.
Operating steps and pre-registered reference numbers are in
[week5_homework2_guide.md](week5_homework2_guide.md); this file reports what the tool
actually produced.

---

## 1. Data and preprocessing

Uploaded via the "RNA-seq Transcriptomics (Course Demo)" one-click loader, which reads the
same local files the assignment names (verified against `webapp/backend/helpers/demo_data.R`,
whose first-priority path for this loader is `tests/RNAseq_output.csv` /
`tests/RNAseq_mapping.csv`):

- **Assay:** 24,393 genes × 24 samples, raw integer counts.
- **Metadata:** `RNAseq_mapping.csv` — a single `Group` column that encodes two crossed
  factors as combined labels (`DMSO`, `DMSO+LIPUS`, `T4400`, `T4400+LIPUS`, `T3976`,
  `T3976+LIPUS`), 3 compounds (vehicle DMSO, and two test compounds T4400/T3976) × 2
  ultrasound conditions (±LIPUS), n = 4 per cell.
- **Import trimming:** the tool drops genes with zero counts in every sample on import —
  24,393 → 19,150 genes. This is a housekeeping step, not a modelling choice.
- **Filter (Preprocess → Filter):** `MIN MAX COUNT = 10` (a gene's peak count across all 24
  samples must reach 10) and `MIN DETECT RATE = 0.167` (≈ 4/24, detected in at least the
  smallest group's worth of samples). This is this tool's closest equivalent to "≥10 counts
  in ≥4 samples"; it is not identical (max-count and detect-rate are evaluated separately,
  and detect-rate only checks count > 0, not count ≥ 10), so it is reported here rather than
  assumed equivalent.
- **Normalisation:** left untouched (no Log Transform applied) before differential testing.
  DESeq2 fits its own size factors from raw counts; a pre-log-transformed matrix would be the
  wrong input for it, which the tool's own AI-interpretation panel flagged unprompted.

## 2. Differential analysis

Run with the **One-click Run All** pipeline (转录组 RNA-seq / DESeq2), reference group
**DMSO**, thresholds |log2FC| ≥ 1 and padj < 0.05 (BH), for each of the three group pairs the
assignment requires. This tool only supports a **flat pairwise comparison** between two
levels of the single `Group` column — it does not fit a joint
`~ compound + lipus + compound:lipus` model — so the compound×ultrasound **interaction is
not testable in this tool**, and no interaction result is claimed here (see §4 for the
independent factorial check).

| Contrast | Genes tested | Significant | Up | Down |
|---|---:|---:|---:|---:|
| **T4400 vs DMSO** | 14,240 | **237** | 170 | 67 |
| T3976 vs DMSO | 14,235 | **1** | 0 | 1 |
| DMSO+LIPUS vs DMSO ("LIPUS within DMSO") | 14,189 | **0** | 0 | 0 |

Only T4400 produces a transcriptional response. T3976's one hit (*Krt19*, log2FC −1.22,
padj 0.040) sits right at the significance boundary and is not evidence of a broad effect.
Ultrasound alone (DMSO+LIPUS vs DMSO) changes nothing detectable.

### Cross-check against an independent factorial model

An independent PyDESeq2 fit of the full factorial design
(`~ compound + lipus + compound:lipus`, script
[`code/hw2_emp_rnaseq_reference.py`](code/hw2_emp_rnaseq_reference.py), numbers in
[week5_homework2_guide.md](week5_homework2_guide.md)) gives **145** significant genes for the
T4400 main effect (82 up, 63 down), 0 for T3976, and 0 for both compound×LIPUS interaction
terms. The web tool's pairwise T4400 count (237) is higher than the joint-model count (145),
in the same direction with the same top genes (below) — the discrepancy is attributable to
method, not disagreement about the biology:

- The web tool compares only the 8 samples in the two groups being contrasted, with no batch
  or ultrasound covariate and (as far as the reported columns show) no apeglm-style shrinkage
  of the fold-change estimates. Smaller n and no shrinkage means noisier log2FC estimates,
  which lets more genes cross a fixed |log2FC| ≥ 1 line even under the same FDR control.
- The joint model borrows dispersion information across all 24 samples and shrinks effect
  sizes, which is more conservative for exactly this reason.

Six of the eight named genes in the reference's "strongest T4400 responses" table appear in
the web tool's result with near-identical fold changes:

| Gene | Web tool log2FC | Reference log2FC | Direction |
|---|---:|---:|---|
| *Mmp13* | +5.23 | +5.25 | Up |
| *Ccl3* | +5.12 | +5.19 | Up |
| *Clec4e* | +4.87 | +4.85 | Up |
| *Kng1* | +4.40 | +4.41 | Up |
| *Lcn2* | +3.51 | +3.52 | Up |
| *Ucma* | −2.55 | −2.54 | Down |

(*Snrpf* and *Mycn*, significant only in the joint model, did not reach padj < 0.05 in the
8-sample pairwise test — consistent with lower power from fewer replicates, not a
contradiction.) A large disagreement here would mean one of the two analyses had the design
wrong; a same-direction, same-magnitude, higher-count-from-less-power pattern like this one
does not, so both are reported rather than one being discarded.

## 3. Functional enrichment (T4400 vs DMSO only — the only contrast with a gene list to enrich)

KEGG and GO enrichment (clusterProfiler, mouse) on the 237-gene T4400 list converge on the
same biology from two independent databases:

**KEGG** — IL-17 signalling pathway (13 genes, padj 5.2×10⁻⁹), cytokine–cytokine receptor
interaction (17 genes, padj 8.6×10⁻⁷), TNF signalling pathway (11 genes, padj 2.4×10⁻⁶),
viral protein interaction with cytokine/cytokine receptor (11 genes, padj 4.8×10⁻⁷), and
rheumatoid arthritis (9 genes, padj 1.4×10⁻⁵).

**GO (biological process)** — acute-phase response (11 genes, padj 1.9×10⁻⁹), acute
inflammatory response (13 genes, padj 5.0×10⁻⁸), antimicrobial humoral response (15 genes,
padj 1.4×10⁻⁷), regulation of body fluid levels (22 genes, padj 3.9×10⁻⁹).

## 4. Interpretation

T4400 drives a coherent innate-immune/inflammatory transcriptional programme — cytokine and
chemokine signalling (IL-17, TNF, cytokine–cytokine receptor interaction), acute-phase and
acute inflammatory response genes — while the two ultrasound-related contrasts (T3976 and
LIPUS alone) show essentially no transcriptional response at this sample size. Among the
individual genes, the strongest inductions are matrix-degrading and myeloid/inflammatory
(*Mmp13*, a collagenase; *Ccl3* and *Clec4e*, myeloid chemokine/pattern-recognition genes;
*Lcn2* and *Kng1*, acute-phase genes), alongside suppression of *Ucma*, a cartilage-matrix
gene expressed by differentiated chondrocytes. Read together with the enrichment results
(rheumatoid arthritis pathway genes, acute inflammatory response), the defensible statement
is that T4400 induces an inflammatory, matrix-remodelling transcriptional shift, consistent
with the chondrocyte catabolic phenotype used to model joint tissue degeneration — not a
claim about T4400's specific molecular target, which this dataset does not identify.

**What this dataset does not show:**

1. **No interaction test.** This tool does not test compound×ultrasound interaction; the
   independent factorial cross-check in §2 found neither interaction term significant
   (smallest padj 0.96 for T4400×LIPUS, 0.999 for T3976×LIPUS), so there is no evidence here
   that ultrasound changes what either compound does — but that claim rests on the external
   check, not on anything computed inside EasyMultiProfiler-Web.
2. **T3976's null result is a power statement, not a biological one.** One borderline gene
   out of 14,235 tested, with n = 4 per group, means this design cannot resolve a modest
   T3976 effect if one exists — it does not mean T3976 is inert.
3. **Association, not mechanism.** Nothing here shows T4400 acts directly on chondrocytes
   rather than another cell type in the tissue, and there is no dose or time course.
4. **Count vs joint-model disagreement (§2)** is expected and explained by design
   (pairwise vs factorial, no covariate control, no shrinkage), not a sign that either result
   is wrong.

## 5. Submission record

- **Synced to GitHub:** T4400 vs DMSO, commit
  [`6186eaa`](https://github.com/santianwan/Bioinformatics_homework_zhouhanyi/commit/6186eaa2153b0ed96ce76fcf59f2a4ed26958cda),
  landed at `EMP2026/Week_05/transcriptomics/weekly/runs/2026-09-24T18-02-32-860Z-dys473/`.
  Sync captures one active analysis state at a time; this run is the headline result (the
  only contrast with a nonzero gene list to enrich and interpret).
- **T3976 vs DMSO and DMSO+LIPUS vs DMSO** were run and downloaded as local result bundles
  (One-click Run All ZIPs) rather than separately synced, since Sync only carries the
  currently active comparison and these two are documented in full in §2 above (1 gene and 0
  genes respectively).
- All three contrasts' DESeq2 output tables, summaries and volcano plots (including the
  synced T4400 run, kept here too for convenience) are checked into
  [`outputs/hw2_webtool_runs/`](outputs/hw2_webtool_runs/) — the numbers in §2–§3 above are
  read directly from these files, not retyped.
