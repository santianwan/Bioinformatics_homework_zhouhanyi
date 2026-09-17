"""Week 4 figures. Run from Week4/code/ ; writes into ../figures/."""
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Arc
import numpy as np, pathlib

FIG = pathlib.Path(__file__).resolve().parent.parent / "figures"
FIG.mkdir(exist_ok=True)

# one palette, threaded across all four figures (§4.1)
C_Q      = "#1b3a5c"   # the biological question / decision nodes
C_DNA    = "#2b7bba"   # DNA-level assays
C_CHROM  = "#5aa469"   # chromatin-level assays
C_RNA    = "#d98c1f"   # RNA / expression
C_FUNC   = "#a3364a"   # functional perturbation (alarm hue, §4.5)
C_GREY   = "#888888"
C_DROP   = "#b9bcc0"


def box(ax, x, y, w, h, text, fc, tc="white", fs=7, style="round,pad=0.012", lw=0, ls="-",
        head=None):
    """head = optional first line drawn in real bold (mathtext would mangle hyphens)."""
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=style, linewidth=lw,
                                facecolor=fc, edgecolor=tc if lw else "none",
                                linestyle=ls, mutation_aspect=1))
    if head is None:
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
                fontsize=fs, color=tc, linespacing=1.35, zorder=5)
    else:
        nlines = text.count("\n") + 2
        pad = h / (nlines + 0.6)
        ax.text(x + w / 2, y + h - pad * 0.75, head, ha="center", va="center",
                fontsize=fs + .4, fontweight="bold", color=tc, zorder=5)
        ax.text(x + w / 2, y + (h - pad * 1.5) / 2, text, ha="center", va="center",
                fontsize=fs, color=tc, linespacing=1.35, zorder=5)


def arrow(ax, p0, p1, color=C_GREY, rad=0.0, lw=1.0, ls="-"):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=9,
                                 linewidth=lw, color=color, linestyle=ls,
                                 shrinkA=1, shrinkB=1,
                                 connectionstyle=f"arc3,rad={rad}", zorder=3))


def blank(ax):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")


# ---------------------------------------------------------------- Q1 ---------
def fig_q1():
    fig, ax = plt.subplots(figsize=(8.6, 5.6)); blank(ax)
    ax.set_title("A mechanism question is answered by ordering assays, not by running all of them",
                 fontsize=8.5, loc="left", pad=8)

    box(ax, .12, .900, .76, .082,
        "Gene $X$ is up-regulated in disease vs control (RNA-seq).\n"
        "Which layer causes it — sequence, accessibility, TF/histone state, methylation, or 3D contact?",
        C_Q, fs=7)

    y1 = .690
    t1 = [(.045, "ATAC-seq", "open chromatin at the\n$X$ promoter and enhancers", C_CHROM),
          (.365, "WGS / WES", "sequence variants in $X$\nand its cis-regulatory region", C_DNA),
          (.685, "WGBS / EM-seq", "CpG methylation at the\n$X$ promoter and CGI", C_DNA)]
    for x, name, sub, c in t1:
        box(ax, x, y1, .27, .110, f"measures: {sub}", c, fs=6.3, head=name)
        arrow(ax, (.5, .900), (x + .135, y1 + .110), C_GREY)
    ax.text(.045, y1 - .030, "run first, in parallel: three cheap, hypothesis-free readouts of the candidate locus",
            fontsize=6, color=C_GREY, style="italic")

    box(ax, .31, .520, .38, .060, "Which layer moved with expression?", C_Q, fs=7)
    for x, _, _, _ in t1:
        arrow(ax, (x + .135, y1), (.5, .580), C_GREY, rad=0.12 if x < .4 else -0.12)

    y2 = .310
    t2 = [(.030, "CUT&Tag / ChIP-seq", "TF or histone-mark\noccupancy", C_CHROM, "if accessibility changed"),
          (.275, "CAGE / RAMPAGE", "which TSS is\nactually fired", C_RNA, "if the promoter is suspect"),
          (.520, "Hi-C / Micro-C", "enhancer-promoter\ncontact frequency", C_CHROM, "if the driver is distal"),
          (.765, "MPRA", "whether a sequence\ncan drive transcription", C_FUNC, "if a variant is suspect")]
    for x, name, sub, c, cond in t2:
        box(ax, x, y2, .205, .110, f"measures:\n{sub}", c, fs=6.0, head=name)
        arrow(ax, (.5, .520), (x + .1025, y2 + .110), C_GREY, rad=0.18 if x < .5 else -0.18)
        ax.text(x + .1025, y2 - .028, cond, fontsize=5.6, color=C_GREY, ha="center", style="italic")

    box(ax, .16, .112, .68, .080,
        "the only step that tests causality: change in $X$ expression\n"
        "when the candidate element is silenced or removed",
        C_FUNC, fs=6.6, head="CRISPRi  /  enhancer deletion  /  base editing")
    for x, _, _, _, _ in t2:
        arrow(ax, (x + .1025, y2), (.5, .192), C_FUNC, rad=0.14 if x < .5 else -0.14, lw=.9)

    box(ax, .16, .012, .68, .056,
        "Mechanism assigned to one layer, with a measured effect size on Gene $X$", C_Q, fs=7)
    arrow(ax, (.5, .112), (.5, .068), C_Q, lw=1.2)

    ax.plot([.045, .955], [.222, .222], color=C_GREY, lw=.7, ls=(0, (3, 3)))
    ax.text(.955, .232, "above this line: correlation only", fontsize=6, color=C_GREY,
            ha="right", va="bottom", style="italic")

    fig.savefig(FIG / "Q1_workflow.png", dpi=300, bbox_inches="tight", facecolor="white")
    return fig


# ---------------------------------------------------------------- Q2 ---------
def fig_q2():
    fig = plt.figure(figsize=(9.0, 6.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.22, 1.0], wspace=.10)
    axL = fig.add_subplot(gs[0]); blank(axL)
    axR = fig.add_subplot(gs[1]); blank(axR)

    axL.set_title("Short-read WGS workflow: every step names its output format",
                  fontsize=8.2, loc="left", pad=8)

    steps = [
        ("1  FASTQ QC",             "FastQC + MultiQC; fastp / Trim Galore", "trimmed FASTQ.gz", C_DNA),
        ("2  Reference genome",     "GRCh38 full analysis set, ALT-aware", "FASTA + .fai + .dict", C_DNA),
        ("3  Alignment",            "BWA-MEM2, with an -R read-group string", "unsorted BAM", C_CHROM),
        ("4  Mapped-read process",  "samtools sort; MarkDuplicates; BQSR", "sorted CRAM + index", C_CHROM),
        ("5  Variant calling",      "GATK HaplotypeCaller, or DeepVariant", "gVCF, then VCF", C_RNA),
        ("6  Annotation",           "VEP + gnomAD + ClinVar, cache pinned", "annotated VCF", C_RNA),
        ("7  Visualization",        "IGV pileups; MultiQC; coverage plots", "PNG / HTML", C_Q),
        ("8  Interpretation",       "ACMG evidence ladder; phenotype gene list", "ranked candidates", C_FUNC),
    ]
    h, gap = .088, .030
    y = .880
    for name, tool, fmt, c in steps:
        box(axL, .075, y, .64, h, tool, c, fs=6.2, head=name)
        axL.text(.735, y + h / 2, fmt, fontsize=5.8, color="#333333", va="center")
        if y > .08:
            arrow(axL, (.395, y), (.395, y - gap), C_GREY, lw=.9)
        y -= (h + gap)

    axL.add_patch(FancyArrowPatch((.075, .880 + h / 2), (.075, .880 - (h + gap) * 3 + h / 2),
                                  arrowstyle="-|>", mutation_scale=8, lw=.9, color=C_FUNC,
                                  connectionstyle="arc3,rad=0.45", zorder=3))
    axL.text(.008, .690, "re-QC after\ntrimming;\nre-align if\nQC fails", fontsize=5.6,
             color=C_FUNC, ha="left", va="center", linespacing=1.35)

    axR.set_title("S01 demo FASTQ, measured (120 PE reads $\\times$ 80 bp)",
                  fontsize=8.2, loc="left", pad=8)
    obs = [
        ("Adapter content", "AGATCGGAAGAGC in 18/120 reads (15.0%) of\nboth mates, always starting at base 68",
         "FAIL. Trim. The fixed start reveals a\nsynthetic insert, not real read-through.", C_FUNC),
        ("Per-base quality", "mean Q 36.0 (cycles 1-10) falls to 27.5 (71-80);\n18/120 reads have median tail Q = 12",
         "WARN. A 15% read subset crashes,\nnot the whole lane.", C_RNA),
        ("Per-sequence GC", "R1 mean 42.8% GC, 12/120 reads (10.0%) > 65%;\nR2 mean 41.3% GC, none > 65%",
         "WARN. The high-GC spike is R1-only, so\nit cannot be a real fragment contaminant.", C_RNA),
        ("Duplication", "R1: 86 unique of 120, top template x18 (15.0%);\nR2: 120 unique of 120, top template x1",
         "WARN. R1-only duplicate: paired\nMarkDuplicates would not flag it.", C_FUNC),
        ("Per-base N", "0 N bases in 9,600 bases per mate",
         "PASS. No basecalling failure.", C_CHROM),
        ("Length distribution", "a single 80 bp bin in both mates",
         "PASS. Uniform, as expected before trimming.", C_CHROM),
    ]
    yy = .880
    hh = .134
    for name, val, call, c in obs:
        axR.add_patch(Rectangle((.010, yy - hh + .018), .011, hh - .018, facecolor=c, edgecolor="none"))
        axR.text(.040, yy, name, fontsize=6.6, fontweight="bold", color="#1a1a1a", va="top")
        axR.text(.040, yy - .030, val, fontsize=5.7, color="#333333", va="top", linespacing=1.45)
        axR.text(.040, yy - .080, call, fontsize=5.7, color=c, va="top", linespacing=1.45,
                 style="italic")
        yy -= hh
    axR.text(.010, .035, "Measured by Week4/code/q2_fastq_metrics.py directly from the provided\n"
                         "demo_fastq/*.fastq.gz, not copied from fastqc_snapshot.tsv.",
             fontsize=5.6, color=C_GREY, va="bottom", linespacing=1.5)

    fig.savefig(FIG / "Q2_workflow.png", dpi=300, bbox_inches="tight", facecolor="white")
    return fig


# ---------------------------------------------------------------- Q3 ---------
def fig_q3():
    fig = plt.figure(figsize=(7.2, 5.6))
    gs = fig.add_gridspec(6, 1, height_ratios=[1, 1, 1, 1.25, 1, 1.15], hspace=.30)
    x = np.linspace(0, 100, 1000)
    ELEM, PROM = 22.0, 74.0          # candidate element, Gene Y promoter (kb)

    def peak(c, w, a):
        return a * np.exp(-((x - c) ** 2) / (2 * w ** 2))

    tracks = [
        ("ATAC-seq\naccessibility", C_CHROM, peak(ELEM, 1.6, 1.0) + peak(PROM, 1.3, .78) + .03),
        ("H3K27ac\nCUT&Tag", C_CHROM, peak(ELEM - 1.6, 2.2, .82) + peak(ELEM + 1.8, 2.2, .90)
         + peak(PROM, 1.8, .55) + .03),
    ]
    axes = []
    for i, (lab, c, y) in enumerate(tracks):
        ax = fig.add_subplot(gs[i]); axes.append(ax)
        ax.fill_between(x, 0, y, color=c, alpha=.85, linewidth=0)
        ax.set_ylim(0, 1.18); ax.set_xlim(0, 100)
        ax.set_yticks([]); ax.set_xticks([])
        for s in ax.spines.values(): s.set_visible(False)
        ax.text(-.012, .5, lab, transform=ax.transAxes, ha="right", va="center",
                fontsize=6.4, linespacing=1.25)

    # methylation
    ax = fig.add_subplot(gs[2]); axes.append(ax)
    rng = np.random.default_rng(4)
    cpg = np.sort(rng.uniform(0, 100, 150))
    meth = np.clip(.82 - .78 * np.exp(-((cpg - ELEM) ** 2) / (2 * 2.2 ** 2))
                   - .74 * np.exp(-((cpg - PROM) ** 2) / (2 * 1.8 ** 2))
                   + rng.normal(0, .05, cpg.size), 0, 1)
    ax.scatter(cpg, meth, s=5, color=C_DNA, linewidth=0)
    ax.set_ylim(-.05, 1.08); ax.set_xlim(0, 100); ax.set_xticks([])
    ax.set_yticks([0, 1]); ax.set_yticklabels(["0", "1"], fontsize=5.6)
    for s in ["top", "right", "bottom"]: ax.spines[s].set_visible(False)
    ax.text(-.012, .5, "WGBS\nCpG methylation", transform=ax.transAxes, ha="right",
            va="center", fontsize=6.4, linespacing=1.25)

    # Hi-C contact arc
    ax = fig.add_subplot(gs[3]); axes.append(ax)
    ax.add_patch(Arc(((ELEM + PROM) / 2, 0), PROM - ELEM, 1.7, theta1=0, theta2=180,
                     color=C_CHROM, lw=2.4))
    ax.add_patch(Arc((ELEM + 14, 0), 28, .72, theta1=0, theta2=180,
                     color=C_DROP, lw=1.2))
    ax.text((ELEM + PROM) / 2, .92, "significant loop\n(element ↔ Gene $Y$ promoter)",
            fontsize=5.8, ha="center", va="bottom", color=C_CHROM, linespacing=1.25)
    ax.text(ELEM + 14, .05, "weaker contact with the intervening gene $Z$ promoter",
            fontsize=5.5, ha="center", va="bottom", color=C_GREY)
    ax.set_ylim(0, 1.5); ax.set_xlim(0, 100); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)
    ax.text(-.012, .5, "Micro-C\n3D contact", transform=ax.transAxes, ha="right",
            va="center", fontsize=6.4, linespacing=1.25)

    # gene track
    ax = fig.add_subplot(gs[4]); axes.append(ax)
    ax.add_patch(Rectangle((ELEM - 1.9, .52), 3.8, .3, facecolor=C_FUNC, edgecolor="none"))
    ax.text(ELEM, .93, "candidate element\n(~4 kb, upstream)", fontsize=6, ha="center",
            va="bottom", color=C_FUNC, linespacing=1.25)
    ax.add_patch(Rectangle((PROM, .60), 19, .14, facecolor=C_RNA, edgecolor="none"))
    ax.annotate("", xy=(PROM + 4.5, .90), xytext=(PROM, .90),
                arrowprops=dict(arrowstyle="-|>", color=C_RNA, lw=1.1))
    ax.text(PROM + 9.5, .30, "Gene $Y$", fontsize=6.6, ha="center", va="top", color=C_RNA)
    ax.add_patch(Rectangle((ELEM + 20, .60), 9, .14, facecolor=C_DROP, edgecolor="none"))
    ax.text(ELEM + 24.5, .30, "Gene $Z$ (alternative target)", fontsize=5.6,
            ha="center", va="top", color=C_GREY)
    ax.set_ylim(0, 1.25); ax.set_xlim(0, 100); ax.set_yticks([])
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xticklabels(["0", "25", "50", "75", "100 kb"], fontsize=5.6)
    for s in ["top", "right", "left"]: ax.spines[s].set_visible(False)
    ax.text(-.012, .5, "Genes", transform=ax.transAxes, ha="right", va="center",
            fontsize=6.4)

    # the logic chain
    ax = fig.add_subplot(gs[5]); blank(ax)
    chain = [("accessible\n(ATAC peak)", C_CHROM), ("active mark\n(H3K27ac)", C_CHROM),
             ("hypo-\nmethylated", C_DNA), ("contacts\nGene $Y$", C_CHROM),
             ("Gene $Y$\nexpressed", C_RNA), ("CRISPRi →\n$Y$ falls?", C_FUNC)]
    w = .142
    for i, (t, c) in enumerate(chain):
        xx = .012 + i * (w + .026)
        box(ax, xx, .40, w, .46, t, c, fs=5.9)
        if i < len(chain) - 1:
            arrow(ax, (xx + w, .63), (xx + w + .026, .63), C_GREY, lw=1.0)
    ax.text(.012, .20, "correlative evidence — consistent with, but not proof of, an enhancer",
            fontsize=5.8, color=C_GREY, style="italic")
    ax.text(.845, .20, "causal test", fontsize=5.8, color=C_FUNC, style="italic",
            ha="center")
    ax.plot([.828, .828], [.34, .92], color=C_GREY, lw=.6, ls=(0, (3, 3)))

    axes[0].set_title("Five correlative layers converge on one element; only perturbation closes the chain",
                      fontsize=8, loc="left", pad=8)
    for a in axes:
        a.axvspan(ELEM - 2.6, ELEM + 2.6, color=C_FUNC, alpha=.07, linewidth=0, zorder=0)
        a.axvspan(PROM - 1.8, PROM + 1.8, color=C_RNA, alpha=.07, linewidth=0, zorder=0)

    fig.savefig(FIG / "Q3_locus_chain.png", dpi=300, bbox_inches="tight",
                facecolor="white")
    return fig


# ---------------------------------------------------------------- Q4 ---------
def fig_q4():
    fig = plt.figure(figsize=(9.2, 5.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.05], wspace=.34)
    ax = fig.add_subplot(gs[0])
    axR = fig.add_subplot(gs[1]); blank(axR)

    # counts read off ../results/q4_filter_trace.tsv (produced by q4_filter.R)
    stages = [
        ("All called variants", 12, None),
        ("FILTER = PASS", 10, "2 dropped: LowQual ($\\it{MSH2}$), FAIL ($\\it{HLA}$-$\\it{A}$)"),
        ("DP $\\geq$ 20 and GQ $\\geq$ 30", 9, "1 dropped: $\\it{MECP2}$, DP=5, GQ=20"),
        ("AF $\\leq$ 0.001", 6, "3 dropped: $\\it{F5}$ 0.42, $\\it{ATM}$ 0.18, chr4 0.35"),
        ("Alters protein or transcript", 4, "2 dropped: $\\it{CFTR}$ synonymous, chr8 intergenic"),
        ("Gene symbol assigned", 4, "none dropped here"),
        ("Not ClinVar (likely) benign", 3, "1 dropped: $\\it{LDLR}$ likely benign"),
    ]
    labs = [s[0] for s in stages]; vals = [s[1] for s in stages]
    ypos = np.arange(len(stages))[::-1]
    for yv, v, (_, _, note) in zip(ypos, vals, stages):
        ax.barh(yv, v, height=.42, color=C_DNA if v > 3 else C_FUNC, linewidth=0)
        ax.text(v + .25, yv, str(v), va="center", fontsize=6.8, color="#1a1a1a")
        if note:
            ax.text(.15, yv - .34, note, fontsize=5.5, color=C_GREY, va="top")
    ax.set_yticks(ypos); ax.set_yticklabels(labs, fontsize=6.4)
    ax.set_xlim(0, 13.6); ax.set_xlabel("variants surviving")
    ax.set_xticks([0, 3, 6, 9, 12])
    ax.set_title("Every drop is recorded, because every filter is a claim",
                 fontsize=8.2, loc="left", pad=8)
    set_frame(ax, "open"); ax.margins(y=.10)

    axR.set_title("The three survivors, after independent annotation checks",
                  fontsize=8.2, loc="left", pad=8)
    cands = [
        ("1", "chr17:7673803 G>A  $\\it{TP53}$",
         "splice acceptor; ClinVar Pathogenic; DP 80, GQ 99, AF 1e-5",
         "Coordinate confirmed inside $\\it{TP53}$ on GRCh38 (Ensembl)\n"
         "and REF=G matches hg38 (UCSC). But VCV000012345\n"
         "resolves to a $\\it{TNFRSF1A}$ record, so the Pathogenic\n"
         "label itself is unverifiable.", C_FUNC),
        ("2", "chr13:32316461 C>T  $\\it{BRCA2}$",
         "missense; ClinVar uncertain significance; DP 60, GQ 90, AF 1e-4",
         "Position is inside $\\it{BRCA2}$ (Ensembl), but the hg38\n"
         "reference base is A, not C, and VCV000067890\n"
         "resolves to $\\it{SCN5A}$. Stays a VUS.", C_RNA),
        ("3", "chr12:25398284 C>A  $\\it{KRAS}$",
         "missense; ClinVar conflicting; DP 58, GQ 91, AF 1.5e-4",
         "A GRCh37 coordinate: GRCh38 $\\it{KRAS}$ spans\n"
         "12:25,205,246-25,326,473. Also a somatic hotspot,\n"
         "not germline biology. De-prioritised.", C_GREY),
    ]
    yy = .900; hh = .268
    for rank, title, qual, verif, c in cands:
        axR.add_patch(Rectangle((.010, yy - hh + .045), .011, hh - .055, facecolor=c, edgecolor="none"))
        axR.text(.040, yy, f"{rank}.  {title}", fontsize=6.8, fontweight="bold",
                 color="#1a1a1a", va="top")
        axR.text(.040, yy - .034, qual, fontsize=5.7, color="#333333", va="top")
        axR.text(.040, yy - .066, verif, fontsize=5.7, color=c, va="top", linespacing=1.5)
        yy -= hh
    axR.text(.010, .000, "Checked against Ensembl REST (gene at coordinate, both builds), the UCSC\n"
                         "hg38/hg19 reference base, and NCBI ClinVar esummary for every stated VCV ID.\n"
                         "Full table: results/q4_annotation_verification.csv",
             fontsize=5.6, color=C_GREY, va="bottom", linespacing=1.5)

    fig.savefig(FIG / "Q4_prioritization.png", dpi=300, bbox_inches="tight", facecolor="white")
    return fig


if __name__ == "__main__":
    apply_figure_style(frame="open", sizes=(8, 6.5, 5.8))
    figs = {"Q1": fig_q1(), "Q2": fig_q2(), "Q3": fig_q3(), "Q4": fig_q4()}
    print("wrote:", sorted(p.name for p in FIG.glob("*.png")))
