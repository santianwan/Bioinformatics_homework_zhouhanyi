"""
Week 5 — PyDESeq2 cross-run of the DESeq2 analysis.

WHY THIS EXISTS. The graded artefact is `week5_deseq2_analysis.R`. It could not be
executed in the analysis sandbox: bioconda ships no win-64 build of DESeq2, and
Bioconductor's own Windows binaries install but refuse to load, because the sandbox
will not dyn.load a compiled .dll from a writable directory. Rather than hand over a
script whose numbers nobody had ever seen, the same design is fitted here with
PyDESeq2 0.5.4 — the reference Python reimplementation of the DESeq2 method, including
apeglm shrinkage — so that every figure in the report is real output.

The R script remains the submission. Where the two disagree, the R run wins.

Run from Week5/code/ :  python week5_pydeseq2_crossrun.py
Outputs -> ../outputs/ and ../figures/  (all prefixed pydeseq2_)
"""
import json
import pathlib

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE.parent / "data"
OUT = HERE.parent / "outputs"
FIG = HERE.parent / "figures"
OUT.mkdir(exist_ok=True, parents=True)
FIG.mkdir(exist_ok=True, parents=True)

PADJ_CUT, LFC_CUT = 0.05, 1.0

# The sandbox cannot create the named pipes joblib's default "loky" backend needs.
# PyDESeq2 always opens its joblib context with `inner_max_num_threads=1`, which
# joblib only accepts from a process-based backend, so plain "threading" and
# "sequential" are both rejected. Register a sequential backend that declares
# support for that argument: the work then runs in-process, one task at a time.
# This changes scheduling only -- every fitted value is identical.
from joblib import register_parallel_backend
from joblib._parallel_backends import SequentialBackend


class _SandboxSequentialBackend(SequentialBackend):
    supports_inner_max_num_threads = True


register_parallel_backend("sandbox_sequential", _SandboxSequentialBackend)
INFERENCE = DefaultInference(n_cpus=1, backend="sandbox_sequential")

# ------------------------------------------------------------------ 1. Import
counts = pd.read_csv(DATA / "Week5_Homework_Count_Matrix.csv", index_col=0)
meta = pd.read_csv(DATA / "Week5_Homework_Sample_Metadata.csv", index_col=0)

# ------------------------------------------------------- 2. Same validation
assert counts.shape[1] == meta.shape[0], "sample count mismatch"
assert list(counts.columns) == list(meta.index), "sample IDs differ or out of order"
assert not counts.columns.duplicated().any(), "duplicated sample IDs"
assert not counts.index.duplicated().any(), "duplicated gene IDs"
assert not counts.isna().any().any(), "missing values"
assert (counts.values >= 0).all(), "negative counts"
assert np.array_equal(counts.values, np.round(counts.values)), "non-integer counts"
assert (counts.values == 0).any(), "no zeros - matrix may not be raw counts"

lib = counts.sum(axis=0)
assert (lib.max() - lib.min()) / lib.mean() > 1e-6, "column sums constant - normalised?"

meta["condition"] = pd.Categorical(meta["condition"], categories=["control", "treated"])
meta["batch"] = pd.Categorical(meta["batch"], categories=["A", "B", "C"])

print("=== Design ===")
print(pd.crosstab(meta["batch"], meta["condition"]))
print(f"\nlibrary sizes: {lib.min():,} - {lib.max():,}  ({lib.max()/lib.min():.2f}x)")
print(f"genes detected per sample: {(counts > 0).sum(axis=0).min()} - "
      f"{(counts > 0).sum(axis=0).max()} of {counts.shape[0]}\n")

# ----------------------------------------------------------- 3/4. Filter, then build
# Filtering happens on the count table, before the object is constructed: slicing a
# DeseqDataSet returns a plain AnnData and silently loses the DESeq2 methods.
# Rule: keep genes with >=10 counts in >=3 samples.
n_before = counts.shape[0]
keep = (counts >= 10).sum(axis=1) >= 3
counts_f = counts.loc[keep]
print(f"=== Filtering: >=10 counts in >=3 samples ===\n"
      f"genes before: {n_before}\ngenes after:  {counts_f.shape[0]} "
      f"({100 * (1 - counts_f.shape[0] / n_before):.1f}% removed)\n")

# PyDESeq2 takes samples as rows. `condition` is an ordered Categorical with
# control first, so the formula treats control as the reference level and the
# fitted coefficient is condition[T.treated].
dds = DeseqDataSet(
    counts=counts_f.T,
    metadata=meta,
    design="~batch + condition",
    refit_cooks=True,
    inference=INFERENCE,
    quiet=True,
)

dm = dds.obsm["design_matrix"]
assert np.linalg.matrix_rank(dm.values) == dm.shape[1], "design matrix rank deficient"
print(f"model matrix: {dm.shape[0]} samples x {dm.shape[1]} coefficients, "
      f"rank {np.linalg.matrix_rank(dm.values)} (full rank)")
print("coefficients:", list(dm.columns), "\n")

# ------------------------------------------------------------------- 5. Fit
dds.deseq2()

# ------------------------------------------------- 6/7. Contrast, then shrink
stat = DeseqStats(dds, contrast=["condition", "treated", "control"],
                  alpha=PADJ_CUT, inference=INFERENCE, quiet=True)
stat.summary()
res_raw = stat.results_df.copy()

coeff = [c for c in dm.columns if c.startswith("condition")]
assert len(coeff) == 1, f"expected one condition coefficient, got {coeff}"
coeff = coeff[0]
print(f"shrinking coefficient: {coeff}")

stat.lfc_shrink(coeff=coeff)
res = stat.results_df.copy()

# apeglm must not move the adjusted p values; if it did, the wrong coefficient
# was shrunk.
assert np.allclose(res["padj"].fillna(-1), res_raw["padj"].fillna(-1), atol=1e-10), \
    "shrinkage altered adjusted p values - wrong coefficient?"
print(f"shrinkage: median |log2FC| {res_raw['log2FoldChange'].abs().median():.3f} "
      f"-> {res['log2FoldChange'].abs().median():.3f}")

# Direction check by hand, independent of the library's labelling.
norm = pd.DataFrame(dds.layers["normed_counts"], index=dds.obs_names, columns=dds.var_names)
top = res["padj"].idxmin()
mt = norm.loc[meta["condition"] == "treated", top].mean()
mc = norm.loc[meta["condition"] == "control", top].mean()
manual = np.log2(mt / mc)
print(f"direction check on {top}: treated {mt:.1f} vs control {mc:.1f}, "
      f"manual log2FC {manual:+.2f}, reported {res.loc[top, 'log2FoldChange']:+.2f}")
assert np.sign(manual) == np.sign(res.loc[top, "log2FoldChange"]), "sign disagrees"

# --------------------------------------------------------------- 8. Result table
res = res.reset_index().rename(columns={"index": "gene_id"})
res["significant"] = res["padj"].notna() & (res["padj"] < PADJ_CUT) & \
                     (res["log2FoldChange"].abs() >= LFC_CUT)
res["direction"] = np.where(res["significant"] & (res["log2FoldChange"] > 0), "Up in treated",
                   np.where(res["significant"] & (res["log2FoldChange"] < 0), "Down in treated",
                            "Not significant"))
res = res.sort_values("padj", na_position="last")
res.to_csv(OUT / "pydeseq2_results.csv", index=False)

n_sig = int(res["significant"].sum())
n_up = int((res["direction"] == "Up in treated").sum())
n_dn = int((res["direction"] == "Down in treated").sum())
print(f"\n=== padj < {PADJ_CUT} AND |log2FC| >= {LFC_CUT} ===")
print(res["direction"].value_counts().to_string())
print(f"padj < {PADJ_CUT} regardless of effect size: "
      f"{int((res['padj'] < PADJ_CUT).sum())}")
print(f"padj = NA: {int(res['padj'].isna().sum())}\n")

# ------------------------------------------------------------------- 9. VST + PCA
dds.vst(use_design=True)
vst = pd.DataFrame(dds.layers["vst_counts"], index=dds.obs_names, columns=dds.var_names)

# DESeq2's plotPCA uses the 500 most variable genes.
top_var = vst.var(axis=0).nlargest(min(500, vst.shape[1])).index
X = vst[top_var].values
Xc = X - X.mean(axis=0)
U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
pcs = U * S
pct = 100 * S**2 / (S**2).sum()

pca_df = pd.DataFrame({"PC1": pcs[:, 0], "PC2": pcs[:, 1],
                       "condition": meta["condition"].astype(str).values,
                       "batch": meta["batch"].astype(str).values},
                      index=dds.obs_names)
pca_df.to_csv(OUT / "pydeseq2_pca_coords.csv")
print(f"=== PCA ===\nPC1 {pct[0]:.0f}% variance, PC2 {pct[1]:.0f}%")
print(pca_df.groupby("condition", observed=True)[["PC1", "PC2"]].mean().to_string())
print(pca_df.groupby("batch", observed=True)[["PC1", "PC2"]].mean().to_string())

from scipy import stats as sstats
g = [pca_df.loc[pca_df.condition == c, "PC1"].values for c in ("control", "treated")]
p_cond_pc1 = sstats.f_oneway(*g).pvalue
gb1 = [pca_df.loc[pca_df.batch == b, "PC1"].values for b in ("A", "B", "C")]
gb2 = [pca_df.loc[pca_df.batch == b, "PC2"].values for b in ("A", "B", "C")]
p_batch_pc1 = sstats.f_oneway(*gb1).pvalue
p_batch_pc2 = sstats.f_oneway(*gb2).pvalue
print(f"\nPC1 ~ condition: p = {p_cond_pc1:.3g}")
print(f"PC1 ~ batch:     p = {p_batch_pc1:.3g}")
print(f"PC2 ~ batch:     p = {p_batch_pc2:.3g}\n")

summary = dict(
    engine=f"PyDESeq2 {__import__('pydeseq2').__version__}",
    genes_input=int(n_before), genes_tested=int(dds.n_vars), samples=int(dds.n_obs),
    design="~batch + condition", coefficient=coeff,
    padj_cut=PADJ_CUT, lfc_cut=LFC_CUT,
    n_significant=n_sig, n_up=n_up, n_down=n_dn,
    n_padj_only=int((res["padj"] < PADJ_CUT).sum()),
    n_padj_na=int(res["padj"].isna().sum()),
    lib_min=int(lib.min()), lib_max=int(lib.max()),
    pc1_var=round(float(pct[0])), pc2_var=round(float(pct[1])),
    p_pc1_condition=float(p_cond_pc1), p_pc1_batch=float(p_batch_pc1),
    p_pc2_batch=float(p_batch_pc2), top_gene=str(top),
)
(OUT / "pydeseq2_run_summary.json").write_text(json.dumps(summary, indent=2))
print("wrote outputs/pydeseq2_results.csv, pydeseq2_pca_coords.csv, "
      "pydeseq2_run_summary.json")
