"""
Week 4, Question 4 — verify the annotation columns of variants_q4.tsv against
primary resources instead of trusting them.

Three independent checks per variant:
  1. Ensembl REST     which gene overlaps the stated coordinate, in GRCh38 AND GRCh37
  2. UCSC REST        the actual reference base at the stated position, hg38 and hg19
  3. NCBI ClinVar     what the stated VCV accession actually describes

The point is not "is this table real" (README_data.md already says it is synthetic).
The point is: if the GENE / REF / CLINVAR_ID columns had been trusted, what would
have been wrong? Answers land in ../results/q4_annotation_verification.csv.

Needs network access. Ensembl returns intermittent 500s under load, hence the retry.

Run:  python q4_verify_annotations.py
"""
import csv
import json
import pathlib
import time
import urllib.error
import urllib.parse
import urllib.request

HERE = pathlib.Path(__file__).resolve().parent
VARIANTS = HERE.parent / "data" / "variants_q4.tsv"
OUT = HERE.parent / "results"
OUT.mkdir(exist_ok=True)

ENSEMBL_38 = "https://rest.ensembl.org"
ENSEMBL_37 = "https://grch37.rest.ensembl.org"
UCSC = "https://api.genome.ucsc.edu"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

PAUSE = 0.4          # be polite to shared public APIs
RETRIES = 6
UNAVAILABLE = "unavailable (API error)"

# Ensembl and UCSC both return intermittent 5xx under load, and a failed lookup must
# never be scored as a mismatch: "the API did not answer" and "the annotation is
# wrong" are different findings. Successful responses are cached on disk so a re-run
# retries only what actually failed.
CACHE_PATH = HERE.parent / "results" / ".api_cache.json"
try:
    _CACHE = json.loads(CACHE_PATH.read_text(encoding="utf-8"))
except (OSError, ValueError):
    _CACHE = {}


def get_json(url):
    if url in _CACHE:
        return _CACHE[url]
    last = None
    for attempt in range(RETRIES):
        try:
            req = urllib.request.Request(
                url, headers={"Accept": "application/json",
                              "Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=60) as fh:
                payload = json.load(fh)
            time.sleep(PAUSE)
            _CACHE[url] = payload
            CACHE_PATH.write_text(json.dumps(_CACHE), encoding="utf-8")
            return payload
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as exc:
            last = exc
            time.sleep(2.0 * (attempt + 1))
    return {"__error__": f"{type(last).__name__}: {last}"}


def genes_at(server, chrom, pos):
    c = chrom.replace("chr", "")
    j = get_json(f"{server}/overlap/region/human/{c}:{pos}-{pos}"
                 f"?feature=gene;content-type=application/json")
    if isinstance(j, dict):
        return UNAVAILABLE
    named = [g.get("external_name") or f"[unnamed {g.get('biotype')}]" for g in j]
    return "; ".join(named) if named else "no gene at this position"


def gene_span(server, symbol):
    j = get_json(f"{server}/lookup/symbol/homo_sapiens/{symbol}"
                 f"?content-type=application/json")
    if "start" not in j:
        return UNAVAILABLE
    return f"{j['seq_region_name']}:{j['start']}-{j['end']}"


def ref_base(db, chrom, pos):
    """UCSC getData/sequence is 0-based half-open, so one base is [pos-1, pos)."""
    j = get_json(f"{UCSC}/getData/sequence?genome={db};chrom={chrom};"
                 f"start={pos - 1};end={pos}")
    if not isinstance(j, dict) or "dna" not in j:
        return UNAVAILABLE
    return j["dna"].upper()


def clinvar_records(vcv_ids):
    """VCV000012345 -> uid 12345. esummary returns the record's title, which names
    the transcript and gene the accession actually belongs to."""
    uids = {i: i.replace("VCV", "").lstrip("0") for i in vcv_ids}
    if not uids:
        return {}
    j = get_json(f"{EUTILS}/esummary.fcgi?db=clinvar&retmode=json&id="
                 + ",".join(uids.values()))
    result = j.get("result", {}) if isinstance(j, dict) else {}
    out = {}
    for vcv, uid in uids.items():
        entry = result.get(uid)
        out[vcv] = (entry.get("title") if isinstance(entry, dict)
                    else "NOT RETURNED BY NCBI")
    return out


def main():
    with open(VARIANTS, encoding="utf-8") as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    rows = list(csv.DictReader(lines, delimiter="\t"))
    print(f"{len(rows)} variants to verify\n")

    cv = clinvar_records([r["CLINVAR_ID"] for r in rows
                          if r["CLINVAR_ID"] not in (".", "")])

    spans = {}
    for sym in {r["GENE"] for r in rows if r["GENE"] != "."}:
        spans[sym] = (gene_span(ENSEMBL_38, sym), gene_span(ENSEMBL_37, sym))

    out_rows = []
    for r in rows:
        chrom, pos = r["CHROM"], int(r["POS"])
        g38 = genes_at(ENSEMBL_38, chrom, pos)
        g37 = genes_at(ENSEMBL_37, chrom, pos)
        b38, b19 = ref_base("hg38", chrom, pos), ref_base("hg19", chrom, pos)
        s38, s37 = spans.get(r["GENE"], ("n/a", "n/a"))

        g38_ok = g38 != UNAVAILABLE
        g37_ok = g37 != UNAVAILABLE
        gene_ok_38 = g38_ok and r["GENE"] != "." and r["GENE"] in g38
        gene_ok_37 = g37_ok and r["GENE"] != "." and r["GENE"] in g37
        if not (g38_ok and g37_ok):
            build = "not resolved (API error)"
        elif gene_ok_38 and b38 == r["REF"]:
            build = "GRCh38"
        elif gene_ok_37 and not gene_ok_38:
            build = "GRCh37"
        else:
            build = "indeterminate / synthetic"

        cvrec = cv.get(r["CLINVAR_ID"], "n/a (no ID given)")
        out_rows.append({
            "CHROM": chrom, "POS": pos, "stated_GENE": r["GENE"], "stated_REF": r["REF"],
            "Ensembl_GRCh38_gene_at_POS": g38,
            "Ensembl_GRCh37_gene_at_POS": g37,
            "Ensembl_GRCh38_span_of_stated_gene": s38,
            "Ensembl_GRCh37_span_of_stated_gene": s37,
            "UCSC_hg38_ref_base": b38, "UCSC_hg19_ref_base": b19,
            "REF_matches_GRCh38": ("not resolved" if b38 == UNAVAILABLE
                                   else ("yes" if b38 == r["REF"] else "no")),
            "gene_matches_GRCh38": ("not resolved" if not g38_ok
                                    else ("yes" if gene_ok_38 else "no")),
            "coordinate_build_inferred": build,
            "stated_CLINVAR_ID": r["CLINVAR_ID"],
            "stated_CLINVAR_SIG": r["CLINVAR_SIG"],
            "NCBI_ClinVar_record_for_that_ID": cvrec,
            "ID_matches_stated_gene": (
                "n/a" if r["CLINVAR_ID"] in (".", "")
                else ("not resolved" if cvrec == "NOT RETURNED BY NCBI"
                      else ("yes" if r["GENE"] and f"({r['GENE']})" in str(cvrec)
                            else "no"))),
        })
        print(f"{r['GENE'] or '.':8s} {chrom}:{pos:<10} GRCh38 gene={g38[:34]:34s} "
              f"REF={r['REF']} hg38={b38} -> {build}")

    with open(OUT / "q4_annotation_verification.csv", "w", newline="",
              encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out_rows[0]))
        w.writeheader()
        w.writerows(out_rows)

    def tally(col, want):
        resolved = [r for r in out_rows if r[col] not in ("not resolved", "n/a")]
        return sum(1 for r in resolved if r[col] == want), len(resolved)

    n_ref, d_ref = tally("REF_matches_GRCh38", "yes")
    n_cv, d_cv = tally("ID_matches_stated_gene", "yes")
    n_37 = sum(1 for r in out_rows if r["coordinate_build_inferred"] == "GRCh37")
    n_unres = sum(1 for r in out_rows
                  if r["coordinate_build_inferred"] == "not resolved (API error)")
    print(f"\nREF matches GRCh38:            {n_ref}/{d_ref} resolved "
          f"({len(out_rows) - d_ref} lookup(s) unavailable)")
    print(f"ClinVar ID matches its gene:  {n_cv}/{d_cv} resolved")
    print(f"Coordinates that are GRCh37:  {n_37}/{len(out_rows)}"
          f"  ({n_unres} not resolved)")
    if n_unres or d_ref < len(out_rows):
        print("NOTE: unresolved lookups are reported as such, never as mismatches. "
              "Re-run to retry only the failed queries (successful ones are cached).")
    print(f"\nwrote {OUT / 'q4_annotation_verification.csv'}")


if __name__ == "__main__":
    main()
