"""
Week 5 — figures and the post-hoc truth benchmark, from the PyDESeq2 cross-run.

Produces:
  ../figures/week5_pca.png               PCA of the VST-transformed samples
  ../figures/week5_de_plot.png           volcano, treated vs control
  ../figures/week5_truth_benchmark.png   recovered vs simulated effect sizes
  ../outputs/week5_truth_benchmark.csv
  ../outputs/week5_benchmark_summary.json

On the instructor key: `Week5_Homework_Gene_Annotation_Instructor_Key.csv` carries
`truth_log2FC_for_instructor`, the simulation's real effect sizes. It was NOT used to
choose any threshold, filter, model term or gene list — week5_deseq2_analysis.R and
week5_pydeseq2_crossrun.py never open it. It enters only here, after the analysis was
fixed, to ask how well the workflow recovered a signal whose answer is known.

Run from Week5/code/ :  python week5_figures_and_benchmark.py
"""
import json
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats as sstats

HERE = pathlib.Path(__file__).resolve().parent
DATA = HERE.parent / "data"
OUT = HERE.parent / "outputs"
FIG = HERE.parent / "figures"
FIG.mkdir(exist_ok=True, parents=True)

PADJ_CUT, LFC_CUT = 0.05, 1.0
C_CTRL, C_TREAT = "#2F6DB3", "#C0392B"
C_MISS, C_FP, C_GREY = "#D98C1F", "#6C3483", "#9AA0A6"

try:
    apply_figure_style(frame="open", sizes=(9, 7.5, 6.5))   # noqa: F821
except NameError:
    plt.rcParams.update({"font.size": 9, "figure.dpi": 300})

res = pd.read_csv(OUT / "pydeseq2_results.csv")
pca = pd.read_csv(OUT / "pydeseq2_pca_coords.csv", index_col=0)
run = json.loads((OUT / "pydeseq2_run_summary.json").read_text())

# ----------------------------------------------------------------------- PCA
fig, ax = plt.subplots(figsize=(6.4, 4.6))
markers = {"A": "o", "B": "s", "C": "^"}
for cond, colour in (("control", C_CTRL), ("treated", C_TREAT)):
    for b, mk in markers.items():
        sub = pca[(pca.condition == cond) & (pca.batch == b)]
        if sub.empty:
            continue
        ax.scatter(sub.PC1, sub.PC2, c=colour, marker=mk, s=58,
                   edgecolor="white", linewidth=.6, zorder=3,
                   label=f"{cond}, batch {b}")
for name, row in pca.iterrows():
    ax.annotate(name, (row.PC1, row.PC2), textcoords="offset points",
                xytext=(6, 4), fontsize=6, color="#333333")

ax.axvline(0, color=C_GREY, lw=.5, ls=(0, (3, 3)), zorder=1)
ax.set_xlabel(f"PC1: {run['pc1_var']}% variance")
ax.set_ylabel(f"PC2: {run['pc2_var']}% variance")
ax.set_title("Condition separates the samples on PC1; batch does not",
             loc="left", fontsize=9)

# One legend keyed by colour (condition) and one by marker (batch), built by hand
# so the reader is not made to parse six combined entries.
from matplotlib.lines import Line2D
h_cond = [Line2D([], [], marker="o", ls="", color=C_CTRL, label="control"),
          Line2D([], [], marker="o", ls="", color=C_TREAT, label="treated")]
h_batch = [Line2D([], [], marker=m, ls="", color=C_GREY, label=f"batch {b}")
           for b, m in markers.items()]
leg1 = ax.legend(handles=h_cond, loc="upper left", frameon=False, fontsize=7,
                 handletextpad=.4, borderpad=.2)
ax.add_artist(leg1)
ax.legend(handles=h_batch, loc="lower left", frameon=False, fontsize=7,
          handletextpad=.4, borderpad=.2)
ax.margins(0.16)
try:
    set_frame(ax, "open")   # noqa: F821
except NameError:
    pass
fig.savefig(FIG / "week5_pca.png", dpi=300, bbox_inches="tight", facecolor="white")
plt.close(fig)

# ------------------------------------------------------------------- Volcano
fig, ax = plt.subplots(figsize=(6.4, 4.8))
d = res.dropna(subset=["padj"]).copy()
d["nlp"] = -np.log10(d["padj"].clip(lower=1e-300))
groups = [("Not significant", C_GREY, .45, 8),
          ("Down in treated", C_CTRL, .85, 16),
          ("Up in treated", C_TREAT, .85, 16)]
for label, colour, alpha, size in groups:
    sub = d[d.direction == label]
    ax.scatter(sub.log2FoldChange, sub.nlp, s=size, c=colour, alpha=alpha,
               linewidth=0, label=f"{label} (n = {len(sub)})", zorder=2)

ax.axhline(-np.log10(PADJ_CUT), color="#555555", lw=.7, ls="--", zorder=1)
ax.axvline(-LFC_CUT, color="#555555", lw=.7, ls="--", zorder=1)
ax.axvline(LFC_CUT, color="#555555", lw=.7, ls="--", zorder=1)

# Only the single strongest hit is labelled. The gene IDs are synthetic
# (Gene0001 ...) and carry no information, so a cloud of them would be noise;
# one label anchors the narrative to a specific point.
hit = d[d.significant].nsmallest(1, "padj").iloc[0]
ax.annotate(hit.gene_id, (hit.log2FoldChange, hit.nlp),
            textcoords="offset points", xytext=(-8, 6), ha="right",
            fontsize=7, color="#222222")

ax.set_xlabel("Shrunken log2 fold change (treated / control)")
ax.set_ylabel("-log10 adjusted p value")
ax.set_title(f"{run['n_significant']} of {run['genes_tested']} tested genes pass both "
             f"thresholds", loc="left", fontsize=9)
ax.text(0.99, 0.02,
        f"dashed: adjusted p = {PADJ_CUT} and |log2FC| = {LFC_CUT:.0f}",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=6.5, color="#555555")
ax.legend(loc="upper left", frameon=False, fontsize=7, handletextpad=.4)
ax.margins(0.06)
try:
    set_frame(ax, "open")   # noqa: F821
except NameError:
    pass
fig.savefig(FIG / "week5_de_plot.png", dpi=300, bbox_inches="tight", facecolor="white")
plt.close(fig)

# --------------------------------------------------------- Truth benchmark
truth = (pd.read_csv(DATA / "Week5_Homework_Gene_Annotation_Instructor_Key.csv")
         .rename(columns={"truth_log2FC_for_instructor": "true_lfc"}))
bench = res.merge(truth[["gene_id", "true_lfc"]], on="gene_id", how="inner")
assert len(bench) == len(res), "join lost genes"
filtered_out = truth[~truth.gene_id.isin(res.gene_id)]

bench["truly_changed"] = bench.true_lfc.abs() >= LFC_CUT
tp = int((bench.truly_changed & bench.significant).sum())
fp = int((~bench.truly_changed & bench.significant).sum())
fn = int((bench.truly_changed & ~bench.significant).sum())
tn = int((~bench.truly_changed & ~bench.significant).sum())

called_true = bench[bench.significant & bench.truly_changed]
sign_ok = int((np.sign(called_true.log2FoldChange) == np.sign(called_true.true_lfc)).sum())

slope, intercept, r_val, _, _ = sstats.linregress(bench.true_lfc, bench.log2FoldChange)
med_err = float((bench.log2FoldChange - bench.true_lfc).abs().median())
rmse = float(np.sqrt(((bench.log2FoldChange - bench.true_lfc) ** 2).mean()))

summary = dict(
    tested=len(bench), filtered_out=len(filtered_out),
    filtered_out_truly_changed=int((filtered_out.true_lfc.abs() >= LFC_CUT).sum()),
    truly_changed=int(bench.truly_changed.sum()), called=int(bench.significant.sum()),
    tp=tp, fp=fp, fn=fn, tn=tn,
    sensitivity=round(tp / (tp + fn), 3) if tp + fn else None,
    precision=round(tp / (tp + fp), 3) if tp + fp else None,
    empirical_fdr=round(fp / (tp + fp), 3) if tp + fp else None,
    specificity=round(tn / (tn + fp), 3) if tn + fp else None,
    sign_correct=f"{sign_ok}/{len(called_true)}",
    pearson_r=round(float(r_val), 3), slope=round(float(slope), 3),
    intercept=round(float(intercept), 3),
    median_abs_error=round(med_err, 3), rmse=round(rmse, 3),
)
(OUT / "week5_benchmark_summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))

bench["class"] = np.select(
    [bench.truly_changed & bench.significant,
     bench.truly_changed & ~bench.significant,
     ~bench.truly_changed & bench.significant],
    ["Correctly called", "Missed (false negative)", "False positive"],
    default="Correctly not called")
bench.to_csv(OUT / "week5_truth_benchmark.csv", index=False)

fig, ax = plt.subplots(figsize=(6.4, 5.0))
lim = [bench.true_lfc.min() - .3, bench.true_lfc.max() + .3]
ax.plot(lim, lim, ls="--", color="#555555", lw=.8, zorder=1)
for label, colour, size, alpha in (("Correctly not called", C_GREY, 8, .40),
                                   ("Missed (false negative)", C_MISS, 20, .90),
                                   ("False positive", C_FP, 20, .90),
                                   ("Correctly called", C_TREAT, 20, .90)):
    sub = bench[bench["class"] == label]
    ax.scatter(sub.true_lfc, sub.log2FoldChange, s=size, c=colour, alpha=alpha,
               linewidth=0, label=f"{label} (n = {len(sub)})", zorder=2)
for v in (-LFC_CUT, LFC_CUT):
    ax.axvline(v, color="#BBBBBB", lw=.6, ls=":", zorder=1)
    ax.axhline(v, color="#BBBBBB", lw=.6, ls=":", zorder=1)

ax.set_xlabel("True log2 fold change (simulation)")
ax.set_ylabel("Shrunken log2 fold change (apeglm)")
ax.set_title(f"Shrinkage recovers direction but compresses magnitude "
             f"(slope {slope:.2f})", loc="left", fontsize=9)
ax.text(0.99, 0.03, f"dashed: y = x;  Pearson r = {r_val:.2f}",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=6.5, color="#555555")
ax.legend(loc="upper left", frameon=False, fontsize=7, handletextpad=.4)
ax.margins(0.05)
try:
    set_frame(ax, "open")   # noqa: F821
except NameError:
    pass
fig.savefig(FIG / "week5_truth_benchmark.png", dpi=300, bbox_inches="tight",
            facecolor="white")
plt.close(fig)

print("\nwrote figures/week5_pca.png, week5_de_plot.png, week5_truth_benchmark.png")
