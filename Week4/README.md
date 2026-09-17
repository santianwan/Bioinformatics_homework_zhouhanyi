# Week 4 — Genomics: Sequencing, Variant Interpretation, and AI-Assisted Analysis

Submission for Week 4 of *Bioinformatics: From Multi-Omics Data to Discovery*.

**Report: [`Week4_Homework_Report.md`](Week4_Homework_Report.md)**

Assignment source: `Week 4/Homework/for_student/` in
[xielab2017/Bioinformatics_SUAT_2026_FALL](https://github.com/xielab2017/Bioinformatics_SUAT_2026_FALL).
The four questions are answered in the structure the assignment requires — reasoning
before AI, AI-assisted workflow, verification, final conclusion — with one figure per
question and the required closing statement.

---

## What was actually computed

Two parts of this submission are measurements rather than prose, and both produced
results that contradict the annotations supplied with the exercise. Those
contradictions are the substance of the verification sections.

**Q2 — the demo FASTQ was re-measured, not read off the snapshot.**
`code/q2_fastq_metrics.py` recomputes six FastQC modules from the 120 read pairs in
`data/demo_fastq/`. The planted traps reproduce (15.0% adapter, a 15% read subset with
median tail Q 12, a 10% high-GC shoulder, a template at 18 copies), but two of them are
present in **R1 only**: the high-GC shoulder and the duplicated template. A real
fragment-level contaminant or PCR duplicate appears in both mates, so neither artefact
is coherent as a sequencing failure — and the R1-only duplicate is invisible to
`MarkDuplicates`, which keys on the aligned positions of the pair.

**Q4 — the variant table's annotation columns were checked against primary records.**
`code/q4_verify_annotations.py` queries Ensembl REST (both assemblies), the UCSC REST
API, and NCBI ClinVar. Three findings:

| Check | Result |
|---|---|
| Do the stated ClinVar accessions describe the stated variants? | **No — 0 of 10.** All ten VCV IDs are real, retrievable records, and every one belongs to a different gene (VCV000012345 → *TNFRSF1A*, not *TP53*; VCV000000888 → *ARSB*, not *KRAS*). The `CLINVAR_SIG` column is therefore unverifiable. |
| Are the coordinates all in one genome build? | **No.** `chr7:117199644`, `chr19:11200200` and `chr12:25398284` are GRCh37 positions; in GRCh38 they fall in *ST7*, *DOCK6* and outside *KRAS* respectively. Annotating against GRCh38 without checking would have assigned three variants to the wrong gene, silently. |
| Do the REF alleles match GRCh38? | **3 of 12.** A real VCF with this property would fail `bcftools norm --check-ref`. |

The filter itself (`code/q4_filter.R`) scores every variant against every criterion and
records *why* each one was dropped, rather than chaining `filter()` calls — so
`results/q4_filter_trace.tsv` shows all 12 variants with a reason, and the shortlist is
just its `kept == TRUE` subset.

---

## Layout

```
Week4/
├── README.md                         this file
├── Week4_Homework_Report.md          the graded report
├── figures/
│   ├── Q1_workflow.png               assay-ordering decision tree
│   ├── Q2_workflow.png               WGS workflow + measured QC
│   ├── Q3_locus_chain.png            locus tracks + evidence chain
│   └── Q4_prioritization.png         filter cascade + survivors
├── code/
│   ├── q2_fastq_metrics.py           FastQC-style modules from the demo reads
│   ├── q4_filter.R                   completed Q4 filter (from the course starter)
│   ├── q4_verify_annotations.py      Ensembl / UCSC / ClinVar verification
│   └── make_figures.py               regenerates all four figures
├── results/
│   ├── q2_fastq_observed_metrics.tsv measured value vs snapshot claim, per module
│   ├── q2_fastq_measured.json        raw measured QC values
│   ├── q4_shortlist.tsv              the 3 surviving variants, ranked
│   ├── q4_filter_trace.tsv           all 12 variants, with keep/drop reason
│   └── q4_annotation_verification.csv per-variant checks against all three resources
└── data/                             unmodified copies of the provided inputs
    ├── variants_q4.tsv
    ├── sample_manifest.csv
    ├── fastqc_snapshot.tsv
    └── demo_fastq/S01_CTRL_WGS_R{1,2}.fastq.gz
```

`data/` holds copies of the course-provided files so that every script in `code/` runs
standalone from this folder. They are unmodified; all provided data are teaching
synthetics (see the course package's `data/README_data.md`) and contain no real
patient information.

## Reproducing

```bash
cd Week4/code
python q2_fastq_metrics.py          # -> ../results/q2_fastq_measured.json
python q4_verify_annotations.py     # -> ../results/q4_annotation_verification.csv   (needs network)
Rscript q4_filter.R                 # -> ../results/q4_shortlist.tsv, q4_filter_trace.tsv
python make_figures.py              # -> ../figures/*.png
```

Environment of record: Python 3.13 (stdlib only except matplotlib for the figures),
R 4.5.3 with readr 2.2.0 and dplyr 1.2.1. `q4_verify_annotations.py` needs outbound
access to `rest.ensembl.org`, `grch37.rest.ensembl.org`, `api.genome.ucsc.edu` and
`eutils.ncbi.nlm.nih.gov`; Ensembl returns intermittent 500s under load, which the
script retries.

## AI assistance

Used as the assignment requires — to critique designs already drafted, to run
plan-first passes on the Q2 workflow and Q4 filter, and to separate observations from
interpretations in Q3. What it changed and what was rejected is recorded per question
in the report, and summarised in Appendix A3. Every number in the report was computed
by the scripts in `code/` or read from a primary record; none came from model recall.
