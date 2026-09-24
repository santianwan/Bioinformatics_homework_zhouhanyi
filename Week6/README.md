# Week 6 — 16S Microbiome Analysis in EasyMultiProfiler-Web

**Zhou Hanyi · SUAT24000202**

- **Analysis + hypothesis:** [`week6_hypothesis.md`](week6_hypothesis.md)
- **Independent reference script:** [`code/week6_16s_reference.R`](code/week6_16s_reference.R)
  (also embedded inline in the write-up)
- **Cross-check script against the web tool's own output:** [`code/week6_webtool_crosscheck.R`](code/week6_webtool_crosscheck.R)
- **Synced result:** commit [`4e98dea`](https://github.com/santianwan/Bioinformatics_homework_zhouhanyi/commit/4e98dea72b70379cc139779795d0a06be0c58092)

Assignment source: `Week 6/Homework/Homework-Week6.docx` in
[xielab2017/Bioinformatics_SUAT_2026_FALL](https://github.com/xielab2017/Bioinformatics_SUAT_2026_FALL).
The slide is titled "Homework for week 5/6" but this folder covers Week 6's own task:
run the 16S files from EMP's `tests/` folder through EasyMultiProfiler-Web end to end,
Sync the result, and propose a hypothesis backed by named parameters.

## Result in one line

130 stool samples (IBS/UC × before/after × poor/great responders), genus level: on the
raw table, disease is the one factor that significantly separates both alpha diversity
(Kruskal-Wallis p = 0.019) and community composition (PERMANOVA R² = 2.9%, p = 0.001),
and UC's Shannon diversity trends down further after treatment while IBS does not
(p = 0.057 within UC, p = 0.80 within IBS). Re-running the same tests on
EasyMultiProfiler-Web's own exported numbers keeps every direction but pushes both
p-values above 0.05 — its default prevalence/detect-rate filter and genus collapse
remove some of the low-abundance signal along with the noise, which is itself a
finding worth reporting (see §"Checking that against..." in the write-up). No baseline
taxon distinguishes good from poor responders after FDR correction, on either version.

## Required-submission checklist

| Item | Status |
|---|---|
| Full analysis run in EasyMultiProfiler-Web (import → preprocess → diversity/PCoA) | ✅ see §"What was done" in the write-up |
| Submitted via Sync | ✅ commit link in `week6_hypothesis.md` §Submission record |
| Scientific hypothesis, with parameters named | ✅ `week6_hypothesis.md` |
| Independent verification | ✅ `code/week6_16s_reference.R`, cross-checked against the web tool's own numbers |

## Layout

```
Week6/
├── README.md                          this file
├── week6_hypothesis.md                analysis, hypothesis, and both scripts (embedded)
├── code/
│   ├── week6_16s_reference.R          independent vegan-based diversity + PERMANOVA check
│   └── week6_webtool_crosscheck.R     re-runs the same tests on the web tool's own output
├── data/
│   ├── 16S_level-7.csv                unmodified copy of the course input
│   └── 16S_mapping.csv
├── figures/
│   ├── week6_pcoa.png
│   └── week6_alpha_boxplot.png
└── outputs/
    ├── week6_alpha_diversity.csv
    ├── week6_baseline_response_taxa.csv
    └── webtool_run_20260925-024344/   EMP-Web's own alpha/beta tables, plots, summary log
```

HW6's own EasyMultiProfiler-Web sync output (not duplicated here) lands in
`EMP2026/Week_06/microbiome_16s/weekly/runs/` at the repo root.
