#!/usr/bin/env python3
"""Collapse genes into independent loci, at five merge distances.

What it does
    A human geneticist counts loci, not genes. Neighbouring paralogues and
    genes in linkage disequilibrium inflate any gene-level count, and the core
    contains both. This script merges genes on the same chromosome whose
    annotated intervals are closer than a given gap into one locus, and carries
    the layer membership, breadth and endpoint signal of the locus as the
    maximum over its genes.

    Five gaps are used: 0 (overlapping genes only), 100 kb, 250 kb, 500 kb and
    1 Mb. 250 kb is the primary definition; 100 kb is what the paper reports in
    the text alongside the sensitivity table, because it is the conventional
    locus width in the BMD literature.

    The honest finding lives here too: at locus level the *direct* contrast
    between the two layers loses significance on fracture (p ~ 0.20), because
    the site-specific layer collapses to too few loci to contrast. The enrichment
    against background survives at every gap. Both are printed and written out.

Reads
    results/gene_layers.tsv       (stage 02, script 01)

Writes
    results/loci.tsv              one row per gene: its locus id at each gap
    results/locus_table_<gap>kb.tsv   one row per locus, at each gap
    results/locus_summary.tsv     loci, core loci, specific loci, per gap

    The enrichment of Mendelian genes across loci - Figure 1D, the OR against
    merge distance - is computed in 04_analysis/01_enrichment.py, which is the
    first script that also has the truth sides. Nothing here needs them.

Numbers in the paper
    This script feeds the paper. The panel tables it ends up in are
    the panel tables in results/.
    No value is repeated here: a number written into a docstring goes
    stale the first time the analysis is re-run. Every claim the paper
    makes is declared in 06_manuscript/numbers.json with the table and
    column it comes from, and `make verify` recomputes all of them.

Runtime
    Under a minute.

Usage
    python 02_layers/02_locus_clumping.py
    python 02_layers/02_locus_clumping.py --gaps 0 250000
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "00_setup"))
import common  # noqa: E402
import config as cfg  # noqa: E402


def assign_loci(d: pd.DataFrame, gap: int) -> pd.Series:
    """Merge genes on one chromosome whose intervals are within `gap` bp.

    Genes are walked in coordinate order; a new locus starts whenever the
    chromosome changes or the next gene starts more than `gap` past the running
    end of the current locus. The running end is the maximum end seen so far,
    so that one long gene bridges the ones nested inside it.
    """
    order = d.sort_values(["chr", "gene_start"])
    locus_id, current, last_chr, last_end = [], -1, None, None
    for _, row in order.iterrows():
        if row["chr"] != last_chr or row["gene_start"] - last_end > gap:
            current += 1
            last_chr, last_end = row["chr"], row["gene_end"]
        else:
            last_end = max(last_end, row["gene_end"])
        locus_id.append(current)
    return pd.Series(locus_id, index=order.index, name=f"locus_{gap}")


def build_locus_table(d: pd.DataFrame, locus: pd.Series,
                      endpoints: list[str]) -> pd.DataFrame:
    """One row per locus: its size, its layer, its breadth, its endpoint signal."""
    t = d.join(locus)
    key = locus.name
    g = t.groupby(key)
    out = pd.DataFrame({
        "n_genes": g.size(),
        "chr": g["chr"].first(),
        "start": g["gene_start"].min(),
        "end": g["gene_end"].max(),
        "core": g["layer"].apply(lambda s: int((s == "core").any())),
        "specific": g["layer"].apply(lambda s: int((s == "specific").any())),
        "breadth": g["breadth"].max(),
        "peak": g["peak"].max(),
        "genes": g.apply(lambda x: ",".join(sorted(x.index)), include_groups=False),
    })
    for e in endpoints:
        out[e] = g[e].max()
    # A locus that contains a core gene is a core locus, whatever else it holds.
    out.loc[out.core == 1, "specific"] = 0
    out["layer"] = np.where(out.core == 1, "core",
                            np.where(out.specific == 1, "specific", "rest"))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--gaps", type=int, nargs="*", default=cfg.CLUMP_GAPS,
                    help="merge distances in bp (default: 0 100k 250k 500k 1M)")
    args = ap.parse_args()

    d = common.load_layers()
    Z = common.load_gene_matrix()
    endpoints = [e for e in ("frak", "heel", "oto") if e in Z.columns]
    for e in endpoints:
        d[e] = Z[e]

    common.banner("Loci instead of genes")
    print(f"  genes {len(d):,} | core {int((d.layer == 'core').sum())} "
          f"| site-specific {int((d.layer == 'specific').sum())}")
    print(f"  endpoints carried: {', '.join(endpoints) if endpoints else '(none)'}")

    loci_per_gene = pd.DataFrame(index=d.index)
    summary = []
    tables = {}
    print(f"\n{'gap':>8s}{'loci':>8s}{'core loci':>11s}{'specific loci':>15s}"
          f"{'genes/core locus':>19s}")
    for gap in args.gaps:
        locus = assign_loci(d, gap)
        loci_per_gene[locus.name] = locus
        table = build_locus_table(d, locus, endpoints)
        tables[gap] = table
        core = table[table.layer == "core"]
        spec = table[table.layer == "specific"]
        med = core.n_genes.median() if len(core) else np.nan
        print(f"{gap // 1000:6d}kb{len(table):8,d}{len(core):11d}{len(spec):15d}"
              f"{med:19.0f}")
        summary.append({"gap_bp": gap, "n_loci": len(table), "n_core_loci": len(core),
                        "n_specific_loci": len(spec),
                        "median_genes_per_core_locus": float(med) if len(core) else np.nan,
                        "max_genes_per_core_locus": int(core.n_genes.max()) if len(core) else 0})
        common.write_result(table, f"locus_table_{gap // 1000}kb.tsv")

    primary = tables[cfg.CLUMP_PRIMARY_GAP] if cfg.CLUMP_PRIMARY_GAP in tables \
        else tables[args.gaps[-1]]
    print(f"\n=== At {cfg.CLUMP_PRIMARY_GAP // 1000} kb, the primary definition ===")
    core = primary[primary.layer == "core"]
    spec = primary[primary.layer == "specific"]
    rest = primary[primary.layer == "rest"]
    print(f"  loci {len(primary):,} | core {len(core)} | site-specific {len(spec)} "
          f"| rest {len(rest):,}")

    if endpoints:
        print("\n=== Endpoint signal at locus level ===")
        print(f"{'endpoint':12s}{'core':>9s}{'specific':>10s}{'rest':>9s}"
              f"{'core vs specific p':>21s}")
        rows = []
        for e in endpoints:
            p = (stats.mannwhitneyu(core[e].dropna(), spec[e].dropna())[1]
                 if len(core) > 2 and len(spec) > 2 else np.nan)
            print(f"{cfg.ENDPOINTS.get(e, e):12s}{core[e].median():9.2f}"
                  f"{spec[e].median():10.2f}{rest[e].median():9.2f}{p:21.3f}")
            rows.append({"endpoint": e, "core_median_z": round(float(core[e].median()), 3),
                         "specific_median_z": round(float(spec[e].median()), 3),
                         "rest_median_z": round(float(rest[e].median()), 3),
                         "p_core_vs_specific": p})
        common.write_result(pd.DataFrame(rows), "locus_endpoint_contrast.tsv",
                            index=False)
        print("\n  Read this honestly: the layer *contrast* on fracture does not")
        print("  reach significance once genes are collapsed to loci, because the")
        print("  site-specific layer leaves too few loci to contrast against. The")
        print("  enrichment against background does survive; see 04_analysis/01.")

    print("\n=== Dose-response at locus level (breadth = maximum in the locus) ===")
    print(f"{'breadth':>8s}{'loci':>8s}" +
          "".join(f"{cfg.ENDPOINTS.get(e, e)[:10]:>12s}" for e in endpoints))
    for k in range(7):
        s = primary[primary.breadth == k]
        if len(s) < 20:
            continue
        print(f"{k:8d}{len(s):8,d}" +
              "".join(f"{s[e].median():12.2f}" for e in endpoints))

    common.write_result(loci_per_gene, "loci.tsv")
    common.write_result(pd.DataFrame(summary), "locus_summary.tsv", index=False)
    print("\nNow run:  python 02_layers/03_threshold_sweep.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
