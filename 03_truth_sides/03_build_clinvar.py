#!/usr/bin/env python3
"""Build the third truth side: ClinVar, at variant level, with its own control.

What it does
    Streams the ClinVar VCF and counts, per gene, how many variants are
    classified pathogenic or likely pathogenic, how many of those carry a
    skeletal diagnosis in their CLNDN field, and how many are classified
    benign.

    The third count is the point. PanelApp and the HPO are both gene-level
    curations from related pipelines; ClinVar is variant-level and curated
    independently, which makes it real triangulation rather than a third look
    at the same list. And the benign set is a built-in negative control: if
    breadth were picking up "genes people have sequenced a lot", benign
    variants would enrich in the core exactly as pathogenic ones do. They do
    not - OR 1.10, p 0.80 - and that contrast is reported in the paper.

Reads
    $SKELBREADTH_DATA/reference/clinvar.vcf.gz

Writes
    results/truth_clinvar.tsv     gene, n_pathogenic, n_pathogenic_skeletal,
                                  n_benign, and the three binary indicators
    results/clinvar_provenance.txt  release line and record counts

Numbers in the paper
    This script feeds the paper. The panel tables it ends up in are
    the panel tables in results/.
    No value is repeated here: a number written into a docstring goes
    stale the first time the analysis is re-run. Every claim the paper
    makes is declared in 06_manuscript/numbers.json with the table and
    column it comes from, and `make verify` recomputes all of them.

Runtime
    Three to eight minutes, streaming ~4.5 million records. Memory stays flat.

Usage
    python 03_truth_sides/03_build_clinvar.py
"""

from __future__ import annotations

import argparse
import collections
import gzip
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "00_setup"))
import common  # noqa: E402
import config as cfg  # noqa: E402

# CLNSIG values are free-ish text; these match the substrings that matter and
# nothing else. "Conflicting_interpretations" contains neither and is excluded,
# which is intentional: a conflicted variant is not evidence either way.
RE_PATHOGENIC = re.compile(r"CLNSIG=[^;]*(Pathogenic|Likely_pathogenic)", re.I)
RE_BENIGN = re.compile(r"CLNSIG=[^;]*Benign", re.I)
RE_GENE = re.compile(r"GENEINFO=([^;:]+)")
RE_DISEASE = re.compile(r"CLNDN=([^;]+)")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default="truth_clinvar.tsv")
    args = ap.parse_args()

    common.banner("Truth side 3 of 4: ClinVar")

    path = cfg.require(cfg.REFERENCE / "clinvar.vcf.gz", "ClinVar VCF (GRCh37)")
    re_skeletal = re.compile(cfg.CLINVAR_SKELETAL_PATTERN, re.I)

    pathogenic = collections.Counter()
    skeletal = collections.Counter()
    benign = collections.Counter()
    header_lines, n_records = [], 0

    with gzip.open(path, "rt", errors="ignore") as fh:
        for line in fh:
            if line.startswith("#"):
                if line.startswith("##fileDate") or line.startswith("##source"):
                    header_lines.append(line.strip())
                continue
            n_records += 1
            if n_records % 500_000 == 0:
                print(f"    {n_records:,} records", end="\r", flush=True)
            fields = line.split("\t")
            if len(fields) < 8:
                continue
            info = fields[7]
            gene = RE_GENE.search(info)
            if not gene:
                continue
            symbol = gene.group(1)
            if RE_PATHOGENIC.search(info):
                pathogenic[symbol] += 1
                disease = RE_DISEASE.search(info)
                if disease and re_skeletal.search(disease.group(1)):
                    skeletal[symbol] += 1
            elif RE_BENIGN.search(info):
                benign[symbol] += 1

    print(f"    {n_records:,} records read")
    print(f"\n  genes with a pathogenic variant     : {len(pathogenic):,}")
    print(f"  genes with a skeletal diagnosis     : {len(skeletal):,}")
    print(f"  genes with a benign variant         : {len(benign):,}")

    genes = sorted(set(pathogenic) | set(benign))
    d = pd.DataFrame({
        "gene": genes,
        "n_pathogenic": [pathogenic.get(g, 0) for g in genes],
        "n_pathogenic_skeletal": [skeletal.get(g, 0) for g in genes],
        "n_benign": [benign.get(g, 0) for g in genes],
    }).set_index("gene")
    d["has_pathogenic"] = (d.n_pathogenic > 0).astype(int)
    d["has_pathogenic_skeletal"] = (d.n_pathogenic_skeletal > 0).astype(int)
    d["has_benign"] = (d.n_benign > 0).astype(int)

    print("\n=== The ten genes with the most skeletal-diagnosis variants ===")
    top = d.sort_values("n_pathogenic_skeletal", ascending=False).head(10)
    for gene, row in top.iterrows():
        print(f"  {gene:12s}{int(row.n_pathogenic_skeletal):6,d} skeletal "
              f"of {int(row.n_pathogenic):6,d} pathogenic")

    provenance = ["ClinVar truth side", f"records read: {n_records}",
                  f"genes pathogenic: {len(pathogenic)}",
                  f"genes with skeletal diagnosis: {len(skeletal)}",
                  f"genes benign: {len(benign)}",
                  f"skeletal diagnosis pattern: {cfg.CLINVAR_SKELETAL_PATTERN}"]
    provenance += header_lines
    cfg.ensure_results_dir()
    (cfg.RESULTS / "clinvar_provenance.txt").write_text(
        "\n".join(provenance) + "\n", encoding="utf-8")

    common.write_result(d, args.out)
    print("  -> results/clinvar_provenance.txt")
    print("\nNow run:  python 03_truth_sides/04_build_brain_truth.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
