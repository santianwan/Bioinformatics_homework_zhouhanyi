# Week 5 Homework 2 — RNA-seq analysis in EasyMultiProfiler-Web

**Name:** Zhou Hanyi  **Student ID:** SUAT24000202  **Date:** 2026-09-24

> **What I could and could not do.** The upload, the in-browser analysis and the **Sync**
> button all run on your machine at `http://127.0.0.1:8080`. That is unreachable from my
> environment — private/reserved IP addresses are a hard security refusal, and my own
> loopback is a different host — so **you have to click those.** Everything else is done:
> the dataset is identified, the same analysis is fitted independently below so you have
> numbers to check the tool against, and the biological interpretation is drafted.
>
> I did not need your local folder: `tests/` is published in the tool's own repository,
> [xielab2017/EasyMultiProfiler-Web](https://github.com/xielab2017/EasyMultiProfiler-Web)
> (branch `main`), and that is the **only** source of the copies in `data/hw2_emp/`.
> **I have never read your local `EasyMultiProfiler-Web-9.0.4\tests\`** — the host-access
> request for that path was refused — so I have not compared the two and cannot claim they
> match. Upload **your** local copies. If their gene or sample counts differ from the
> figures below, yours are correct and my reference numbers do not apply.

---

## 1. Which files to upload

`tests/` holds four datasets. Only one is RNA-seq:

| File | What it is | Use for this homework? |
|---|---|---|
| **`RNAseq_output.csv`** | 24,394 mouse genes × 24 samples, raw integer counts | ✅ **expression matrix** |
| **`RNAseq_mapping.csv`** | 24 rows: `SampleID`, `Group` | ✅ **sample metadata** |
| `16S_level-7.csv` / `16S_mapping.csv` | microbiome amplicon table | ❌ 16S, not RNA-seq |
| `ChIP/HA_summits_0.05.bed` | ChIP-seq peak summits | ❌ not RNA-seq, and BED not CSV |
| `Clinical*.csv`, `meta-*.csv` | clinical / sample annotation | ❌ not expression data |

So: **`RNAseq_output.csv` + `RNAseq_mapping.csv`.**

## 2. The experimental design

`RNAseq_mapping.csv` gives six groups of four. The single `Group` column actually encodes
**two crossed factors**, and recognising that is the main design decision in this homework:

|  | LIPUS − | LIPUS + |
|---|---|---|
| **DMSO** (vehicle) | 4 | 4 |
| **T4400** | 4 | 4 |
| **T3976** | 4 | 4 |

LIPUS is low-intensity pulsed ultrasound; DMSO is the solvent control; T4400 and T3976 are
two test compounds. A balanced 3 × 2 factorial, n = 4 — so the questions available are the
compound main effect, the ultrasound main effect, and their **interaction** (does
ultrasound change what a compound does?).

If the web tool only offers a flat one-factor group comparison, run the pairwise contrasts
against DMSO and say in your write-up that the interaction was not testable in that mode.
Do not pretend a flat comparison answered the factorial question.

---

## 3. Steps in EasyMultiProfiler-Web

These follow the tool's own README (v9.0.4). Web UI `http://127.0.0.1:8080`, API
`http://127.0.0.1:8000`.

1. **Start it.** `Start-EMP-Web.bat` (or `Run-EMP-Web-Windows.bat`). Closing the browser
   does **not** stop the server — use `Stop-EMP-Web-Windows.bat` when you are finished.
2. **Analyze → Import.** Upload `RNAseq_output.csv` as the expression/count matrix and
   `RNAseq_mapping.csv` as the sample metadata. Confirm 24 samples and 24,394 features
   are recognised, and that the `Group` column is picked up as the grouping variable.
3. **Preprocess.** Apply a low-count filter. I used ≥ 10 counts in ≥ 4 samples — 4 being
   the smallest group size — which keeps 13,842 of 24,394 genes (43.3 % removed). Whatever
   the tool's default is, **write the rule down**; the checklist from Homework 1 asks for it.
4. **Differential analysis.** Reference level **DMSO**. Run at least:
   `T4400 vs DMSO`, `T3976 vs DMSO`, and `LIPUS vs no LIPUS`. If the tool can fit
   `~ compound + lipus + compound:lipus`, also extract the two **interaction** terms — they
   are the only test of whether ultrasound potentiates a compound, and a pairwise
   group comparison does not answer it. Use adjusted *p* < 0.05 with |log2FC| ≥ 1, the same
   thresholds as Homework 1, and report both.
5. **Visualisation.** PCA plus a volcano or MA plot for each contrast.
6. **Export → Sync.**
   - Register / log in: **student ID `SUAT24000202`**, name (required), password ≥ 8 characters.
   - If it asks for a token, paste a GitHub fine-grained PAT with **Contents: Read and write**
     on `santianwan/Bioinformatics_homework_zhouhanyi`, then **绑定仓库**.
   - Course track: **transcriptomics**. Assignment: **Week 5** (weekly, not project).
   - Click **同步到 GitHub**, then open the commit link it returns and confirm.

### Where it will land

The tool writes its own tree, separate from the `Week5/` folder in this submission:

```
EMP2026/Week_05/transcriptomics/weekly/runs/<timestamp>/
    manifest.json      data/    results/    plots/    teaching/
```

Sync is **incremental** — each run creates a new timestamped folder and nothing existing is
deleted. Your repo already has `EMP2026/Week_01/microbiome_16s/weekly/runs/...` from Week 1,
so the mechanism is one you have used before.

---

## 4. Reference results — check the tool against these

Fitted independently with PyDESeq2 (same method as Homework 1), design
`~ compound + lipus + compound:lipus`, filter ≥ 10 counts in ≥ 4 samples, thresholds
adjusted *p* < 0.05 and |log2FC| ≥ 1. Script: [`code/hw2_emp_rnaseq_reference.py`](code/hw2_emp_rnaseq_reference.py).

Exact agreement is not expected — the web tool may filter, normalise or shrink differently.
Large disagreement means one of you has the design wrong, and that is worth finding.

| Contrast | What it asks | Significant | Up | Down | padj < 0.05 alone |
|---|---|---:|---:|---:|---:|
| **T4400 vs DMSO** | compound effect, no ultrasound | **145** | 82 | 63 | 1,610 |
| T3976 vs DMSO | compound effect, no ultrasound | 0 | 0 | 0 | 1 |
| LIPUS within DMSO | ultrasound effect **in the vehicle arm only** | 0 | 0 | 0 | 0 |
| T4400 × LIPUS | does ultrasound change the size of the T4400 effect? | 0 | 0 | 0 | 0 |
| T3976 × LIPUS | does ultrasound change the size of the T3976 effect? | 0 | 0 | 0 | 0 |

Because the model carries an interaction, the `lipus` coefficient is the ultrasound effect
**at the reference compound (DMSO)**, not an overall ultrasound effect — hence the row name.
The last two rows are the actual potentiation test; the smallest adjusted *p* in either is
0.96 (T4400 × LIPUS) and 0.999 (T3976 × LIPUS), so neither is close.

**Strongest T4400 responses**

| Up in T4400 | log2FC | padj | | Down in T4400 | log2FC | padj |
|---|---:|---:|---|---|---:|---:|
| *Mmp13* | +5.25 | 0.026 | | *Snrpf* | −6.93 | 0.033 |
| *Ccl3* | +5.19 | 0.027 | | *Ucma* | −2.54 | 0.009 |
| *Clec4e* | +4.85 | 0.046 | | *1600027J07Rik* | −2.34 | 0.011 |
| *Kng1* | +4.41 | 0.031 | | *Fam198a* | −2.20 | 0.009 |
| *Lcn2* | +3.52 | 0.046 | | *Mycn* | −2.10 | 0.039 |

**QC observations worth reporting**

- **Sequencing depth is bimodal.** Six samples sit near 20.7 M assigned reads and eighteen
  near 40.2 M — almost exactly 2×. This looks alarming, so I tested whether it tracks the
  design: it does not (χ², *p* = 0.75; the six low-depth samples are spread across five of
  the six groups). Median-of-ratios normalisation handles it. **Report the check, not just
  the observation** — a depth split that *had* been confounded with group would have
  invalidated the whole analysis.
- **PC1 carries 80 % of the variance but only 46 % of PC1 is between-group.** Within-group
  scatter on PC1 (SD 3.8–17.1) is comparable to the spread between group means (−8.4 to
  +13.5). The samples do not form tidy clusters, and no sample is more than 1.5 SD from its
  own group mean — the heterogeneity is spread through the cohort, not caused by one
  outlier. **Do not drop any sample.**

---

## 5. Draft interpretation

> Do not paste this verbatim — it is the reasoning, and the marks are for your version of it.

Of the three comparisons, only **T4400** produces a transcriptional response: 145 genes at
adjusted *p* < 0.05 with at least a two-fold change, 82 up and 63 down. **T3976 is
essentially silent** (1 gene at FDR alone, none passing both thresholds), and **ultrasound
changes nothing detectable in the vehicle arm** — no gene passes even the FDR cut for LIPUS
versus no LIPUS within DMSO. The PCA agrees: the only samples displaced along PC1 are
T4400, and LIPUS explains nothing on that axis (*p* = 0.64).

The T4400 gene list is coherent rather than scattered. The strongest inductions —
*Mmp13*, *Ccl3*, *Clec4e*, *Kng1*, *Lcn2* — are a matrix-degradation and innate-inflammatory
set: *Mmp13* is the principal collagenase of cartilage catabolism, *Lcn2* and *Kng1* are
acute-phase, and *Ccl3* and *Clec4e* are myeloid chemokine and pattern-recognition genes.
At the same time *Ucma*, a cartilage-matrix gene expressed by differentiated chondrocytes,
falls ~6-fold. Up-catabolic, down-anabolic, in mouse tissue, with ultrasound as the other
factor: that combination is the transcriptional signature of a **chondrocyte catabolic
shift**, the pattern used to model osteoarthritis-like degeneration. So the defensible
statement is that T4400 drives an inflammatory, matrix-degrading programme while suppressing
cartilage matrix synthesis.

Three limits on that reading, all of which belong in the write-up:

1. **A null is not an absence.** With n = 4 and this much within-group scatter, "T3976 does
   nothing" means *this experiment cannot resolve what T3976 does*. 1,610 genes pass FDR for
   T4400 but only 145 clear the effect-size bar, which is the same point from the other
   direction: the design is underpowered for modest effects.
2. **Neither interaction term is significant.** Fitting
   `~ compound + lipus + compound:lipus` and testing the two `compound:lipus` coefficients
   directly returns no gene at FDR for either compound — the smallest adjusted *p* is 0.96
   for T4400 × LIPUS and 0.999 for T3976 × LIPUS. So there is no evidence here that
   ultrasound potentiates either compound, which was presumably the point of the design.
   Note this is a statement about an interaction test that was actually run: a pairwise
   group comparison cannot support it, and reporting the null honestly is worth more than
   hunting for a sub-threshold trend.
3. **Gene-level association only.** Nothing here shows the compound acts directly on
   chondrocytes rather than through another cell type in the tissue, and no dose or time
   course is available.

---

## 6. Checklist before you sync

- [ ] Both CSVs imported; 24 samples × 24,394 features confirmed
- [ ] Filter rule recorded
- [ ] DMSO set as the reference level
- [ ] Three contrasts run; both thresholds reported
- [ ] PCA + volcano/MA exported
- [ ] Depth-vs-group check mentioned in the interpretation
- [ ] Logged in with student ID `SUAT24000202`; track = transcriptomics; assignment = Week 5
- [ ] Commit link opened and verified after sync
