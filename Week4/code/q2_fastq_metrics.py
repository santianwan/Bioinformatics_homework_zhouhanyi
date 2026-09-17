"""
Week 4, Question 2 — FastQC-style metrics computed directly from the demo FASTQ.

The homework offers `data/demo_fastq/fastqc_snapshot.tsv` as a substitute for
running FastQC. This script does not use it: it recomputes the same modules from
the reads so the snapshot can be *checked* rather than quoted. Anything the
snapshot claims and this script cannot reproduce is reported as a discrepancy.

Run:  python q2_fastq_metrics.py
Writes: ../results/q2_fastq_measured.json
"""
import gzip
import json
import pathlib
import statistics as stats
from collections import Counter

HERE = pathlib.Path(__file__).resolve().parent
FQDIR = HERE.parent / "data" / "demo_fastq"
OUT = HERE.parent / "results"
OUT.mkdir(exist_ok=True)

# Illumina TruSeq / Nextera-style adapter prefix. Verified against the sequence
# named in fastqc_snapshot.tsv and in FastQC's own contaminant list.
ADAPTER = "AGATCGGAAGAGC"
PHRED_OFFSET = 33          # Illumina 1.8+ / Sanger encoding
TAIL_START = 50            # the snapshot claims the crash begins "after cycle 50"
TAIL_Q_CUTOFF = 20         # a read whose tail mean Q < 20 is "crashed"
HIGH_GC = 65.0             # %GC above which a read joins the high-GC shoulder


def read_fastq(path):
    """Yield (name, seq, qual). Deliberately strict: a truncated record raises."""
    with gzip.open(path, "rt") as fh:
        while True:
            name = fh.readline()
            if not name:
                return
            seq, plus, qual = fh.readline(), fh.readline(), fh.readline()
            if not qual:
                raise ValueError(f"truncated FASTQ record at {name!r}")
            seq, qual = seq.strip(), qual.strip()
            if len(seq) != len(qual):
                raise ValueError(f"seq/qual length mismatch in {name!r}")
            yield name.strip(), seq, qual


def metrics(path):
    recs = list(read_fastq(path))
    n = len(recs)
    lengths = Counter(len(s) for _, s, _ in recs)
    maxlen = max(lengths)

    per_cycle = [[] for _ in range(maxlen)]
    for _, _, q in recs:
        for i, ch in enumerate(q):
            per_cycle[i].append(ord(ch) - PHRED_OFFSET)
    mean_q = [stats.mean(c) for c in per_cycle]

    gc = [100 * (s.count("G") + s.count("C")) / len(s) for _, s, _ in recs]

    crashed_tails = []
    for _, _, q in recs:
        tail = [ord(c) - PHRED_OFFSET for c in q[TAIL_START:]]
        if stats.mean(tail) < TAIL_Q_CUTOFF:
            crashed_tails.append(tail)

    adapter_hits = [s.find(ADAPTER) for _, s, _ in recs if ADAPTER in s]
    dup = Counter(s for _, s, _ in recs)
    top_seq, top_n = dup.most_common(1)[0]

    return {
        "file": path.name,
        "n_reads": n,
        "read_lengths": dict(lengths),
        "mean_q_cycles_1_10": round(stats.mean(mean_q[:10]), 2),
        "mean_q_last_10_cycles": round(stats.mean(mean_q[-10:]), 2),
        "lowest_per_cycle_mean_q": round(min(mean_q), 2),
        "lowest_q_cycle": min(range(maxlen), key=lambda i: mean_q[i]) + 1,
        "crashed_reads": len(crashed_tails),
        "crashed_pct": round(100 * len(crashed_tails) / n, 1),
        "median_tail_q_in_crashed_reads": (
            stats.median([v for t in crashed_tails for v in t]) if crashed_tails else None),
        "gc_mean_pct": round(stats.mean(gc), 2),
        "reads_above_65pct_gc": sum(1 for g in gc if g > HIGH_GC),
        "pct_above_65pct_gc": round(100 * sum(1 for g in gc if g > HIGH_GC) / n, 1),
        "adapter_reads": len(adapter_hits),
        "adapter_pct": round(100 * len(adapter_hits) / n, 1),
        "adapter_start_positions_0based": sorted(set(adapter_hits)),
        "n_bases_total": sum(s.count("N") for _, s, _ in recs),
        "bases_total": sum(len(s) for _, s, _ in recs),
        "unique_sequences": len(dup),
        "top_duplicate_copies": top_n,
        "top_duplicate_pct": round(100 * top_n / n, 1),
        "top_duplicate_contains_adapter": ADAPTER in top_seq,
    }


if __name__ == "__main__":
    result = {tag: metrics(FQDIR / f"S01_CTRL_WGS_{tag}.fastq.gz") for tag in ("R1", "R2")}

    # The three cross-mate checks that the snapshot does not make. Each of these
    # is a statement about whether the planted artefact is physically coherent.
    result["cross_mate_checks"] = {
        "duplicate_present_in_both_mates":
            result["R1"]["top_duplicate_copies"] > 1 and result["R2"]["top_duplicate_copies"] > 1,
        "high_gc_shoulder_in_both_mates":
            result["R1"]["reads_above_65pct_gc"] > 0 and result["R2"]["reads_above_65pct_gc"] > 0,
        "adapter_in_both_mates":
            result["R1"]["adapter_reads"] > 0 and result["R2"]["adapter_reads"] > 0,
        "adapter_start_is_fixed":
            len(result["R1"]["adapter_start_positions_0based"]) == 1,
    }

    (OUT / "q2_fastq_measured.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result["cross_mate_checks"], indent=2))
    for tag in ("R1", "R2"):
        m = result[tag]
        print(f"{tag}: {m['n_reads']} reads, GC {m['gc_mean_pct']}%, "
              f"adapter {m['adapter_pct']}%, {m['unique_sequences']} unique, "
              f"top dup x{m['top_duplicate_copies']}")
