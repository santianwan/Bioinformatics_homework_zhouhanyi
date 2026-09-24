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
C = {"DMSO": "#7F8C8D", "T4400": "#2F6DB3", "T3976": "#C0392B"}
MK = {"no": "o", "yes": "^"}

fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.3))
ax = axes[0]
for comp, colour in C.items():
    for lp, mk in MK.items():
        sub = pca[(pca.compound == comp) & (pca.lipus == lp)]
        ax.scatter(sub.PC1, sub.PC2, c=colour, marker=mk, s=52, alpha=.9,
                   edgecolor="white", linewidth=.6, zorder=3)
ax.set_xlabel(f"PC1: {pct[0]:.0f}% variance")
ax.set_ylabel(f"PC2: {pct[1]:.0f}% variance")
ax.set_title("Only T4400 shifts the transcriptome, and within-group scatter is large",
             loc="left", fontsize=9)
from matplotlib.lines import Line2D
h = [Line2D([], [], marker="o", ls="", color=c, label=k) for k, c in C.items()]
h += [Line2D([], [], marker=m, ls="", color="#555555", label=f"LIPUS {k}")
      for k, m in MK.items()]
ax.legend(handles=h, frameon=False, fontsize=7, loc="best", ncol=2)
ax.margins(.12)

ax = axes[1]
names = list(contrasts)
ups = [summary[n]["up"] for n in names]
downs = [-summary[n]["down"] for n in names]
y = np.arange(len(names))
ax.barh(y, ups, height=.55, color="#C0392B", label="Up")
ax.barh(y, downs, height=.55, color="#2F6DB3", label="Down")
for i, (u, d) in enumerate(zip(ups, downs)):
    if u:
        ax.text(u + 8, i, str(u), va="center", fontsize=7, color="#C0392B")
    if d:
        ax.text(d - 8, i, str(-d), va="center", ha="right", fontsize=7, color="#2F6DB3")
    if u == 0 and d == 0:
        # A zero result is a finding here, not missing data: mark it explicitly so
        # the empty slot cannot be read as "not tested".
        ax.plot(0, i, marker="|", color="#333333", ms=9, mew=1.4, zorder=4)
        ax.text(6, i, "0", va="center", fontsize=7, color="#333333")
ax.axvline(0, color="#333333", lw=.8)
ax.set_yticks(y)
ax.set_yticklabels([n.replace("_", " ") for n in names], fontsize=8)
ax.set_xlabel("differentially expressed genes")
ax.set_title(f"Genes passing padj < {PADJ_CUT} and |log2FC| ≥ {LFC_CUT:.0f}",
             loc="left", fontsize=9)
ax.legend(frameon=False, fontsize=7, loc="lower right")
ax.margins(x=.22)
for spine in ("top", "right"):
    for a in axes:
        a.spines[spine].set_visible(False)

fig.tight_layout()
fig.savefig(FIG / "hw2_overview.png", dpi=300, bbox_inches="tight", facecolor="white")
plt.close(fig)

print("\nwrote outputs/hw2_*.csv, hw2_summary.json and figures/hw2_overview.png")
