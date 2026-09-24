# Week 5 — Transcriptomics: RNA-seq and Differential Expression with DESeq2

Submission for Week 5 **Homework 1** of *Bioinformatics: From Multi-Omics Data to Discovery*.

**Zhou Hanyi · SUAT24000202**

- **HW1 interpretation:** [`week5_interpretation.md`](week5_interpretation.md)
- **HW1 AI verification log:** [`week5_AI_verification_log.md`](week5_AI_verification_log.md)
- **HW1 analysis script:** [`code/week5_deseq2_analysis.R`](code/week5_deseq2_analysis.R)
- **HW2 operating guide + independent reference analysis:** [`week5_homework2_guide.md`](week5_homework2_guide.md)
- **HW2 results + interpretation (what the tool actually produced):** [`week5_homework2_results.md`](week5_homework2_results.md)

Assignment source: `Week 5/Homework/for_student/` in
[xielab2017/Bioinformatics_SUAT_2026_FALL](https://github.com/xielab2017/Bioinformatics_SUAT_2026_FALL).

> **Homework 2** (EasyMultiProfiler-Web) runs interactively at `http://127.0.0.1:8080` on the
> student's own machine — the upload, filtering, differential analysis and Sync were done by
> hand there, following [`week5_homework2_guide.md`](week5_homework2_guide.md). The synced
> run (T4400 vs DMSO) lands in `EMP2026/Week_05/transcriptomics/weekly/runs/`, not here;
> [`week5_homework2_results.md`](week5_homework2_results.md) reports the actual numbers from
> that run (and the two comparisons kept as local result bundles) against the independent
> reference analysis.

---

## Result in one line

Treated versus control, `~ batch + condition`, 989 of 1,000 genes tested:
**60 genes significant at adjusted *p* < 0.05 and |log2FC| ≥ 1 — 36 up, 24 down.**
PC1 (24 % of variance) separates condition completely (*p* = 4 × 10⁻¹²); batch
contributes nothing to PC1 (*p* = 0.99).

---

## Read this before running

`week5_deseq2_analysis.R` is the submission. It could not be executed in the environment
where this homework was first drafted (DESeq2 would not install there — see the AI
verification log's original §4 for why), so the reported numbers were first measured with
**PyDESeq2 0.5.4**, the reference Python reimplementation of the DESeq2 method, apeglm
shrinkage included.

**It has since been run on the student's own machine** (R 4.6.1, DESeq2 + apeglm — no
install restriction there). Every number below and in `week5_interpretation.md` is now
confirmed by the R run itself, not inferred from PyDESeq2; see the AI verification log's
2026-09-24 update for the one script fix this run required (`vst()` replaced with
`varianceStabilizingTransformation()`, forced by the filtered gene count falling under
`vst()`'s default subsampling threshold — a fix with no effect on any DE number, only on
how the PCA transform is fit).

### Required-submission checklist

| Required file | Status |
|---|---|
| `week5_deseq2_analysis.R` | ✅ `code/week5_deseq2_analysis.R` |
| `week5_pca.png` | ✅ `figures/` — R render |
| `week5_de_plot.png` | ✅ `figures/` — R render |
| `week5_interpretation.md` | ✅ |
| `week5_AI_verification_log.md` | ✅ |
| `week5_deseq2_results.csv` | ✅ `outputs/week5_deseq2_results.csv` |
| `week5_deseq2_object.rds` | ✅ `outputs/week5_deseq2_object.rds` |
| `session_info.txt` | ✅ `outputs/session_info.txt` |

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
├── week5_homework2_guide.md           HW2 operating steps + independent reference analysis
├── week5_homework2_results.md         HW2 actual results from EasyMultiProfiler-Web + interpretation
├── code/
│   ├── week5_deseq2_analysis.R        THE SUBMISSION - run this locally
│   ├── week5_pydeseq2_crossrun.py     same design fitted with PyDESeq2
│   ├── week5_figures_and_benchmark.py figures + post-hoc truth benchmark
│   ├── week5_truth_benchmark.R        benchmark for the R pipeline
│   └── hw2_emp_rnaseq_reference.py    independent factorial-model cross-check for HW2
├── figures/
│   ├── week5_pca.png
│   ├── week5_de_plot.png
│   └── week5_truth_benchmark.png
├── outputs/
│   ├── pydeseq2_results.csv           full table, 989 genes, non-significant included
│   ├── pydeseq2_pca_coords.csv
│   ├── pydeseq2_run_summary.json
│   ├── week5_truth_benchmark.csv
│   ├── week5_benchmark_summary.json
│   └── hw2_webtool_runs/              EMP-Web's own DESeq2 tables + volcano plots, 3 contrasts
└── data/                              unmodified copies of the provided inputs
```

HW2's own synced results (not duplicated here) land in
`EMP2026/Week_05/transcriptomics/weekly/runs/` at the repo root.

`data/` holds copies of the course files so every script runs standalone from this folder.

## Environment

- R 4.6.1, DESeq2 1.52.0 + apeglm 1.34.0 (`outputs/session_info.txt` has the full `sessionInfo()`).
- Python 3.12, PyDESeq2 0.5.4, pandas, numpy, scipy, matplotlib.

One sandbox-specific note carried in `week5_pydeseq2_crossrun.py`: PyDESeq2 opens its
joblib context with `inner_max_num_threads=1`, which only process-based backends accept,
and the sandbox cannot create the named pipes loky needs. The script registers a
sequential backend that declares support for that argument. This changes scheduling only
— every fitted value is identical.
