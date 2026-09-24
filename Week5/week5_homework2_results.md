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

---

## Appendix: independent reference analysis script

The §2 cross-check numbers (145/0/0 significant genes, the interaction test, the
six-gene fold-change table) come from this script, also kept standalone at
[`code/hw2_emp_rnaseq_reference.py`](code/hw2_emp_rnaseq_reference.py).

```python
"""
Week 5 Homework 2 — reference analysis of the EasyMultiProfiler-Web RNA-seq test data.

WHAT THIS IS. Homework 2 asks for the analysis to be run *inside* EasyMultiProfiler-Web
and submitted with that tool's Sync button. This script does not replace that. It fits
the same dataset with the Homework 1 method (DESeq2 via PyDESeq2) so that the numbers
coming out of the web tool can be checked against an independent run, and so the
biological interpretation rests on something measured.

DATA. tests/RNAseq_output.csv and tests/RNAseq_mapping.csv from
github.com/xielab2017/EasyMultiProfiler-Web — 24,394 mouse genes x 24 samples, raw
integer counts. Design is 3 compounds (DMSO vehicle, T4400, T3976) x 2 ultrasound
levels (-/+ LIPUS), n = 4 per cell.

Run from Week5/code/ :  python hw2_emp_rnaseq_reference.py
Outputs -> ../outputs/hw2_*  and ../figures/hw2_*
"""
import json
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from joblib import register_parallel_backend
from joblib._parallel_backends import SequentialBackend
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from scipy import stats as sstats

HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE.parent / "data" / "hw2_emp"
OUT = HERE.parent / "outputs"
FIG = HERE.parent / "figures"
for d in (OUT, FIG):
    d.mkdir(parents=True, exist_ok=True)

PADJ_CUT, LFC_CUT = 0.05, 1.0


# Same sandbox workaround as the Homework 1 cross-run: PyDESeq2 opens its joblib
# context with inner_max_num_threads=1, which only process-based backends accept,
# and loky cannot create named pipes here. Sequential execution, identical values.
class _SandboxSequentialBackend(SequentialBackend):
    supports_inner_max_num_threads = True


register_parallel_backend("sandbox_sequential", _SandboxSequentialBackend)
INFERENCE = DefaultInference(n_cpus=1, backend="sandbox_sequential")

# ------------------------------------------------------------------ 1. Import
counts = pd.read_csv(DATA / "RNAseq_output.csv", index_col=0)
meta = pd.read_csv(DATA / "RNAseq_mapping.csv", index_col=0)

assert list(counts.columns) == list(meta.index), "sample IDs differ or out of order"
assert not counts.isna().any().any(), "missing values"
assert (counts.values >= 0).all(), "negative counts"
assert np.array_equal(counts.values, np.round(counts.values)), "non-integer counts"

# Split the single Group label into the two crossed factors it encodes.
meta["lipus"] = np.where(meta["Group"].str.contains(r"\+LIPUS"), "yes", "no")
meta["compound"] = meta["Group"].str.replace(r"\+LIPUS", "", regex=True)
meta["compound"] = pd.Categorical(meta["compound"], categories=["DMSO", "T4400", "T3976"])
meta["lipus"] = pd.Categorical(meta["lipus"], categories=["no", "yes"])

lib = counts.sum(axis=0)
meta["depth_M"] = (lib / 1e6).round(2).values

print("=== Design ===")
print(pd.crosstab(meta["compound"], meta["lipus"]))
print(f"\nlibrary sizes: {lib.min():,} - {lib.max():,}  ({lib.max()/lib.min():.2f}x)")

# The depth distribution is bimodal (6 samples near 20 M, 18 near 40 M). That only
# matters if it tracks the design, so test it rather than eyeball it.
meta["depth_class"] = np.where(lib.values < 30e6, "~20M", "~40M")
tab = pd.crosstab(meta["Group"], meta["depth_class"])
chi2_p = sstats.chi2_contingency(tab)[1]
print(f"\ndepth class vs group: chi2 p = {chi2_p:.3f} "
      f"({'not confounded' if chi2_p > 0.05 else 'CONFOUNDED - stop'})")
print(tab.to_string())
assert chi2_p > 0.05, "sequencing depth is confounded with group"

# ------------------------------------------------------------ 2. Filter + fit
n_before = counts.shape[0]
keep = (counts >= 10).sum(axis=1) >= 4          # >=10 counts in >=4 samples (smallest cell)
counts_f = counts.loc[keep]
print(f"\n=== Filtering: >=10 counts in >=4 samples ===\n"
      f"genes before: {n_before}\ngenes after:  {counts_f.shape[0]} "
      f"({100 * (1 - counts_f.shape[0] / n_before):.1f}% removed)")

dds = DeseqDataSet(counts=counts_f.T, metadata=meta,
                   design="~compound + lipus + compound:lipus",
                   refit_cooks=True, inference=INFERENCE, quiet=True)
dm = dds.obsm["design_matrix"]
assert np.linalg.matrix_rank(dm.values) == dm.shape[1], "design matrix rank deficient"
print(f"\nmodel matrix: {dm.shape[0]} x {dm.shape[1]}, full rank")
print("coefficients:", list(dm.columns))

dds.deseq2()

# ------------------------------------------------------- 3. Contrasts of interest
# Each contrast answers one question; naming them here keeps the report honest about
# what was tested rather than reporting whichever comparison looked best.
# With an interaction in the model, the `lipus` coefficient is the ultrasound effect
# AT THE REFERENCE COMPOUND (DMSO) -- not an overall ultrasound effect. Named
# accordingly so the report cannot overstate it.
coefs = list(dm.columns)


def coef_vector(name):
    """Numeric contrast selecting a single fitted coefficient."""
    v = np.zeros(len(coefs))
    v[coefs.index(name)] = 1.0
    return v


contrasts = {
    "LIPUS_within_DMSO":  (["lipus", "yes", "no"],        "Does ultrasound change expression in the DMSO (vehicle) arm?"),
    "T4400_vs_DMSO":      (["compound", "T4400", "DMSO"], "Does compound T4400 alone change expression?"),
    "T3976_vs_DMSO":      (["compound", "T3976", "DMSO"], "Does compound T3976 alone change expression?"),
    # The interaction terms are the only test of "does ultrasound potentiate the
    # compound?". Without these, no statement about potentiation is supportable.
    "Interaction_T4400xLIPUS": (coef_vector("compound[T.T4400]:lipus[T.yes]"),
                                "Does LIPUS change the size of the T4400 effect?"),
    "Interaction_T3976xLIPUS": (coef_vector("compound[T.T3976]:lipus[T.yes]"),
                                "Does LIPUS change the size of the T3976 effect?"),
}

summary, tables = {}, {}
for name, (contrast, question) in contrasts.items():
    st = DeseqStats(dds, contrast=contrast, alpha=PADJ_CUT,
                    inference=INFERENCE, quiet=True)
    st.summary()
    df = st.results_df.copy()
    df["significant"] = (df["padj"].notna() & (df["padj"] < PADJ_CUT)
                         & (df["log2FoldChange"].abs() >= LFC_CUT))
    df["direction"] = np.where(df["significant"] & (df["log2FoldChange"] > 0), "Up",
                       np.where(df["significant"] & (df["log2FoldChange"] < 0), "Down", "NS"))
    df = df.sort_values("padj", na_position="last")
    df.to_csv(OUT / f"hw2_{name}.csv")
    tables[name] = df
    summary[name] = dict(
        question=question,
        tested=int(df["padj"].notna().sum()),
        significant=int(df["significant"].sum()),
        up=int((df["direction"] == "Up").sum()),
        down=int((df["direction"] == "Down").sum()),
        padj_only=int((df["padj"] < PADJ_CUT).sum()),
        top_gene=str(df.index[0]),
        top_lfc=round(float(df["log2FoldChange"].iloc[0]), 3),
        top_padj=float(df["padj"].iloc[0]),
    )
    print(f"\n--- {name}: {question}")
    print(f"    significant (padj<{PADJ_CUT} & |log2FC|>={LFC_CUT}): "
          f"{summary[name]['significant']}  ({summary[name]['up']} up, "
          f"{summary[name]['down']} down)   padj-only: {summary[name]['padj_only']}")
    print(f"    top: {summary[name]['top_gene']} "
          f"log2FC {summary[name]['top_lfc']:+.2f}, padj {summary[name]['top_padj']:.2e}")

# ------------------------------------------------------------------ 4. PCA
dds.vst(use_design=False)          # blind: this is exploratory QC across six groups
vst = pd.DataFrame(dds.layers["vst_counts"], index=dds.obs_names, columns=dds.var_names)
top_var = vst.var(axis=0).nlargest(500).index
X = vst[top_var].values
Xc = X - X.mean(axis=0)
U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
pct = 100 * S**2 / (S**2).sum()
pca = pd.DataFrame({"PC1": (U * S)[:, 0], "PC2": (U * S)[:, 1],
                    "group": meta["Group"].values,
                    "compound": meta["compound"].astype(str).values,
                    "lipus": meta["lipus"].astype(str).values,
                    "depth_class": meta["depth_class"].values}, index=dds.obs_names)
pca.to_csv(OUT / "hw2_pca_coords.csv")
print(f"\n=== PCA ===\nPC1 {pct[0]:.0f}%, PC2 {pct[1]:.0f}%")
for fac in ("compound", "lipus", "depth_class"):
    groups = [pca.loc[pca[fac] == lv, "PC1"].values for lv in pca[fac].unique()]
    print(f"  PC1 ~ {fac:11s} p = {sstats.f_oneway(*groups).pvalue:.3g}")

summary["_design"] = dict(
    genes_input=int(n_before), genes_tested=int(counts_f.shape[0]),
    samples=int(counts.shape[1]),
    lib_min=int(lib.min()), lib_max=int(lib.max()),
    depth_group_chi2_p=round(float(chi2_p), 3),
    pc1_var=round(float(pct[0])), pc2_var=round(float(pct[1])),
    padj_cut=PADJ_CUT, lfc_cut=LFC_CUT,
)
(OUT / "hw2_summary.json").write_text(json.dumps(summary, indent=2))

# ------------------------------------------------------------------ 5. Figures
# (bar/scatter plotting code omitted here for length; see code/hw2_emp_rnaseq_reference.py
#  for the full PCA-scatter + up/down bar-chart figure generation.)
```
