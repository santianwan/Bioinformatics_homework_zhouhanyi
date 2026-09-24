# Week 5 — Transcriptomics: RNA-seq and Differential Expression with DESeq2

Submission for Week 5 **Homework 1** of *Bioinformatics: From Multi-Omics Data to Discovery*.

- **Interpretation:** [`week5_interpretation.md`](week5_interpretation.md)
- **AI verification log:** [`week5_AI_verification_log.md`](week5_AI_verification_log.md)
- **Analysis script:** [`code/week5_deseq2_analysis.R`](code/week5_deseq2_analysis.R)

Assignment source: `Week 5/Homework/for_student/` in
[xielab2017/Bioinformatics_SUAT_2026_FALL](https://github.com/xielab2017/Bioinformatics_SUAT_2026_FALL).

> **Homework 2** (EasyMultiProfiler-Web) is **not** in this folder. It requires uploading
> files through a local web interface and submitting with that tool's own Sync button,
> which has to be done interactively.

---

## Result in one line

Treated versus control, `~ batch + condition`, 989 of 1,000 genes tested:
**60 genes significant at adjusted *p* < 0.05 and |log2FC| ≥ 1 — 36 up, 24 down.**
PC1 (24 % of variance) separates condition completely (*p* = 4 × 10⁻¹²); batch
contributes nothing to PC1 (*p* = 0.99).

---

## Read this before running

`week5_deseq2_analysis.R` is the submission, and **it has not been executed in the
environment where this was prepared.** DESeq2 could not be installed there: bioconda has
no win-64 build, and Bioconductor's Windows binaries install but will not load, because
that sandbox refuses to `dyn.load` a compiled `.dll` from a writable directory.

So that the reported numbers were measured rather than assumed, the identical design was
fitted with **PyDESeq2 0.5.4** — the reference Python reimplementation of the DESeq2
method, apeglm shrinkage included. Every number in the interpretation comes from that
run and is reproducible from `outputs/pydeseq2_results.csv`.

**To produce the R deliverables, run one command on a machine with DESeq2:**

```bash
cd Week5/code
Rscript week5_deseq2_analysis.R
```

That writes `outputs/week5_deseq2_results.csv`, `outputs/week5_deseq2_object.rds`,
`outputs/session_info.txt`, and overwrites `figures/week5_pca.png` and
`figures/week5_de_plot.png` with the R renders. DESeq2 and PyDESeq2 agree closely but not
bit-for-bit; **where they differ, the R output is correct.**

### Required-submission checklist

| Required file | Status |
|---|---|
| `week5_deseq2_analysis.R` | ✅ `code/week5_deseq2_analysis.R` |
| `week5_pca.png` | ✅ `figures/` — PyDESeq2 render; R run overwrites it |
| `week5_de_plot.png` | ✅ `figures/` — PyDESeq2 render; R run overwrites it |
| `week5_interpretation.md` | ✅ |
| `week5_AI_verification_log.md` | ✅ |
| `week5_deseq2_results.csv` | ⬜ **produced by running the R script** (PyDESeq2 equivalent: `outputs/pydeseq2_results.csv`) |
| `week5_deseq2_object.rds` | ⬜ **produced by running the R script** |
| `session_info.txt` | ⬜ **produced by running the R script** |

---

## The part worth reading: how well did the workflow actually work?

The student package ships `Week5_Homework_Gene_Annotation_Instructor_Key.csv`, whose
`truth_log2FC_for_instructor` column holds the effect sizes used to simulate the counts.
It looks like an instructor file.

**It was not used to derive anything** — no threshold, filter, model term or gene list
came from it, and the analysis scripts never open it. It is read only afterwards, to ask
what the exercise otherwise cannot answer: given that the truth is known, how good was
the recovery?

| | |
|---|---|
| Truly changed genes among those tested | 88 |
| Called significant | 60 (58 true positives, 2 false positives) |
| Sensitivity | 0.659 |
| Precision | 0.967 |
| **Empirical FDR** | **0.033** against a nominal 0.05 |
| Direction correct among true positives | **58 / 58** |
| Slope of recovered on true log2FC | 0.843 |

Three conclusions: FDR control is working and is slightly conservative, not
anti-conservative; the direction of effect is never wrong; and **sensitivity, not
specificity, is the binding constraint** — 30 of 88 truly changed genes were missed, all
with true effects near the |log2FC| ≥ 1 cutoff. With six replicates per group this design
is precise and under-powered: what it reports is trustworthy, and it under-reports.

---

## Layout

```
Week5/
├── README.md                          this file
├── week5_interpretation.md            the 142-word interpretation + supporting detail
├── week5_AI_verification_log.md       prompt, what AI got wrong, what I verified
├── code/
│   ├── week5_deseq2_analysis.R        THE SUBMISSION - run this locally
│   ├── week5_pydeseq2_crossrun.py     same design fitted with PyDESeq2
│   ├── week5_figures_and_benchmark.py figures + post-hoc truth benchmark
│   └── week5_truth_benchmark.R        benchmark for the R pipeline
├── figures/
│   ├── week5_pca.png
│   ├── week5_de_plot.png
│   └── week5_truth_benchmark.png
├── outputs/
│   ├── pydeseq2_results.csv           full table, 989 genes, non-significant included
│   ├── pydeseq2_pca_coords.csv
│   ├── pydeseq2_run_summary.json
│   ├── week5_truth_benchmark.csv
│   └── week5_benchmark_summary.json
└── data/                              unmodified copies of the provided inputs
```

`data/` holds copies of the course files so every script runs standalone from this folder.

## Environment

- R 4.5.3, targeting DESeq2 + apeglm (the script's `sessionInfo()` records the actual versions).
- Python 3.12, PyDESeq2 0.5.4, pandas, numpy, scipy, matplotlib.

One sandbox-specific note carried in `week5_pydeseq2_crossrun.py`: PyDESeq2 opens its
joblib context with `inner_max_num_threads=1`, which only process-based backends accept,
and the sandbox cannot create the named pipes loky needs. The script registers a
sequential backend that declares support for that argument. This changes scheduling only
— every fitted value is identical.
