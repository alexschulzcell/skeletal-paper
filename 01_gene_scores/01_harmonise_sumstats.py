#!/usr/bin/env python3
"""Harmonise every GWAS to the two columns MAGMA needs: rsID and p-value.

What it does
    For each trait, reads the archived summary statistics, maps variants to the
    rsIDs of the 1000 Genomes European reference (by chromosome:position where
    the source has no rsID, which is the case for the DXA and GEFOS files), and
    writes a two-column `<trait>.pval` file. Where the source does not state a
    sample size, N is estimated from the standard errors and allele frequencies
    as N = median[ 1 / (2 f (1-f) se^2) ], which is the estimator the GEFOS
    files require and which reproduces their published N to within a percent.

    Three traps this script exists to avoid, all of which cost us a run once:
      - the DXA and GEFOS files carry chr:pos, not rsID, so an unmapped file
        silently yields a MAGMA run over a handful of variants. The script
        aborts below 500,000 mapped variants.
      - a summary file written on Windows leaves \\r on the last field, and
        MAGMA then reads the sample size as "34806\\r" and refuses. Everything
        written here is newline-clean.
      - p-values of exactly 0 (underflow in the source) are dropped rather
        than passed on, where MAGMA would turn them into a silent NaN.

Reads
    $SKELBREADTH_DATA/sumstats/<trait>.tsv.gz
    $SKELBREADTH_DATA/reference/g1000_eur.bim

Writes
    $SKELBREADTH_DATA/magma/<trait>.pval          two columns, SNP and P
    $SKELBREADTH_DATA/magma/sample_sizes.tsv      trait, N, source of N
    results/harmonisation_report.tsv              variants read and mapped

Numbers in the paper
    This script feeds the paper. The panel tables it ends up in are
    the panel tables in results/.
    No value is repeated here: a number written into a docstring goes
    stale the first time the analysis is re-run. Every claim the paper
    makes is declared in 06_manuscript/numbers.json with the table and
    column it comes from, and `make verify` recomputes all of them.

Runtime
    20-60 minutes for all traits; dominated by reading the height file.

Usage
    python 01_gene_scores/01_harmonise_sumstats.py
    python 01_gene_scores/01_harmonise_sumstats.py --traits head fn ls
    python 01_gene_scores/01_harmonise_sumstats.py --force
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "00_setup"))
import config as cfg  # noqa: E402

MIN_MAPPED = 500_000

#: Column names each source family uses, in order of preference.
COLUMN_ALIASES = {
    "rsid": ["variant_id", "rsid", "rs_id", "SNP", "snp", "MarkerName", "rsID"],
    "chr": ["chromosome", "chr", "CHR", "#chrom", "Chr"],
    "pos": ["base_pair_location", "position", "pos", "BP", "bp", "Pos"],
    "p": ["p_value", "p-value", "pval", "P", "p", "P_BOLT_LMM_INF", "P-value"],
    "se": ["standard_error", "se", "SE", "StdErr"],
    "eaf": ["effect_allele_frequency", "eaf", "EAF", "Freq1", "freq"],
    "n": ["n", "N", "sample_size", "n_total"],
}


def pick(columns, kind):
    """First column of the file matching one of the known aliases for `kind`."""
    lower = {c.lower(): c for c in columns}
    for alias in COLUMN_ALIASES[kind]:
        if alias in columns:
            return alias
        if alias.lower() in lower:
            return lower[alias.lower()]
    return None


def load_rsid_map(bim_path: Path) -> pd.Series:
    """chr:pos -> rsID, from the LD reference panel."""
    bim = pd.read_csv(bim_path, sep="\t", header=None, usecols=[0, 1, 3],
                      names=["chr", "snp", "pos"], dtype={"chr": str})
    bim["key"] = bim["chr"] + ":" + bim["pos"].astype(str)
    return bim.drop_duplicates("key").set_index("key")["snp"]


def estimate_n(path: Path, se_col: str, eaf_col: str, limit: int = 400_000) -> int:
    """N from se ~ 1 / sqrt(2 N f (1-f)), over common variants only."""
    se, eaf = [], []
    with gzip.open(path, "rt", errors="ignore") as fh:
        header = fh.readline().strip().replace(",", "\t").split()
        try:
            i_se, i_eaf = header.index(se_col), header.index(eaf_col)
        except ValueError:
            return 0
        for i, line in enumerate(fh):
            if i > limit:
                break
            parts = line.split()
            try:
                f, s = float(parts[i_eaf]), float(parts[i_se])
            except (ValueError, IndexError):
                continue
            if 0.1 < f < 0.9 and s > 0:
                se.append(s)
                eaf.append(f)
    if len(se) < 1000:
        return 0
    se, eaf = np.asarray(se), np.asarray(eaf)
    return int(round(np.median(1.0 / (2 * eaf * (1 - eaf) * se ** 2))))


def harmonise(trait: str, rsid_map: pd.Series, out_dir: Path, force: bool) -> dict:
    src = cfg.SUMSTATS / f"{trait}.tsv.gz"
    dest = out_dir / f"{trait}.pval"
    if dest.exists() and dest.stat().st_size > 0 and not force:
        n_lines = sum(1 for _ in open(dest)) - 1
        print(f"  {trait:12s} already harmonised ({n_lines:,} variants)")
        return {"trait": trait, "read": np.nan, "mapped": n_lines, "status": "cached"}

    if not src.exists():
        print(f"  {trait:12s} SOURCE MISSING ({src.name}) - skipped")
        return {"trait": trait, "read": 0, "mapped": 0, "status": "missing"}

    with gzip.open(src, "rt", errors="ignore") as fh:
        header = fh.readline().strip().replace(",", "\t").split()
    c_rs, c_chr, c_pos, c_p = (pick(header, k) for k in ("rsid", "chr", "pos", "p"))
    if c_p is None:
        print(f"  {trait:12s} NO P-VALUE COLUMN in {header[:8]} - skipped")
        return {"trait": trait, "read": 0, "mapped": 0, "status": "no p column"}

    use = [c for c in (c_rs, c_chr, c_pos, c_p) if c]
    # Sniff the delimiter from the header rather than letting pandas do it with
    # the Python engine, which is several times slower on a ten-million-row file.
    with gzip.open(src, "rt", errors="ignore") as fh:
        first = fh.readline()
    sep = "," if first.count(",") > first.count("\t") else "\t"

    total = mapped = 0
    with open(dest, "w", newline="\n") as out:
        out.write("SNP\tP\n")
        for chunk in pd.read_csv(src, sep=sep, usecols=use,
                                 dtype={c_chr: str} if c_chr else None,
                                 chunksize=2_000_000):
            total += len(chunk)
            if c_rs is not None and chunk[c_rs].astype(str).str.startswith("rs").mean() > 0.5:
                chunk["SNP"] = chunk[c_rs]
            elif c_chr and c_pos:
                key = (chunk[c_chr].astype(str).str.replace("chr", "", regex=False)
                       + ":" + chunk[c_pos].astype("Int64").astype(str))
                chunk["SNP"] = key.map(rsid_map)
            else:
                raise SystemExit(f"{trait}: neither rsID nor chr:pos available")
            chunk = chunk.dropna(subset=["SNP", c_p])
            chunk = chunk[(chunk[c_p] > 0) & (chunk[c_p] <= 1)]
            chunk[["SNP", c_p]].to_csv(out, sep="\t", header=False, index=False,
                                       lineterminator="\n")
            mapped += len(chunk)
            print(f"    {total:>12,} read, {mapped:>12,} mapped", end="\r", flush=True)
    print(" " * 60, end="\r")

    print(f"  {trait:12s} {total:>12,} read, {mapped:>12,} mapped"
          f"{'   *** TOO FEW ***' if mapped < MIN_MAPPED else ''}")
    if mapped < MIN_MAPPED:
        dest.unlink()
        raise SystemExit(
            f"\n{trait}: only {mapped:,} variants mapped to the LD reference "
            f"(threshold {MIN_MAPPED:,}).\nThis is almost always a column-name or "
            f"genome-build mismatch, not a data problem. Header was:\n  {header}\n")
    return {"trait": trait, "read": total, "mapped": mapped, "status": "ok"}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--traits", nargs="*", default=None,
                    help="subset of traits; default is all that are present")
    ap.add_argument("--force", action="store_true",
                    help="re-harmonise traits that already have a .pval file")
    args = ap.parse_args()

    all_traits = (cfg.BMD_SITES + list(cfg.ENDPOINTS)
                  + list(cfg.GEFOS_TRAITS) + cfg.BRAIN_REGIONS)
    traits = args.traits or all_traits

    out_dir = cfg.MAGMA_DIR
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Reading the LD reference to build the chr:pos -> rsID map ...")
    rsid_map = load_rsid_map(cfg.require(cfg.REFERENCE / "g1000_eur.bim",
                                         "1000 Genomes EUR .bim"))
    print(f"  {len(rsid_map):,} reference variants\n")

    print("=== Harmonising ===")
    report = [harmonise(t, rsid_map, out_dir, args.force) for t in traits]

    print("\n=== Sample sizes ===")
    rows = []
    for trait in traits:
        src = cfg.SUMSTATS / f"{trait}.tsv.gz"
        if not src.exists():
            continue
        if trait in cfg.TRAIT_N:
            n, source = cfg.TRAIT_N[trait], "published"
        else:
            with gzip.open(src, "rt", errors="ignore") as fh:
                header = fh.readline().strip().replace(",", "\t").split()
            c_se, c_eaf = pick(header, "se"), pick(header, "eaf")
            n = estimate_n(src, c_se, c_eaf) if (c_se and c_eaf) else 0
            source = "estimated from se and EAF"
        if not n:
            print(f"  {trait:12s} NO SAMPLE SIZE - MAGMA cannot run this trait")
            continue
        print(f"  {trait:12s} N = {n:>10,}   ({source})")
        rows.append({"trait": trait, "N": n, "source": source})

    pd.DataFrame(rows).to_csv(out_dir / "sample_sizes.tsv", sep="\t",
                              index=False, lineterminator="\n")
    cfg.ensure_results_dir()
    pd.DataFrame(report).to_csv(cfg.RESULTS / "harmonisation_report.tsv",
                                sep="\t", index=False)
    print(f"\n  -> {out_dir / 'sample_sizes.tsv'}")
    print(f"  -> results/harmonisation_report.tsv")
    print("\nNow run:  bash 01_gene_scores/02_run_magma.sh")
    return 0


if __name__ == "__main__":
    sys.exit(main())
