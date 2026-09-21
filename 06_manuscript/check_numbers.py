"""Verify every number in the manuscript against the result tables.

The manuscript must never disagree with `results/`. This script extracts each claim
listed in `numbers.json`, recomputes it from the result table it came from, and exits
non-zero on the first mismatch. Run it before every submission and after every change
to the analysis.

    python 06_manuscript/check_numbers.py            # check all claims
    python 06_manuscript/check_numbers.py --verbose  # print every comparison

A claim is defined in numbers.json as:

    {"id": "panelapp_or_breadth6",
     "text": "odds ratio 5.5",              what the manuscript says
     "file": "fig1c_threshold_sweep.csv",   which result table holds the truth
     "where": {"threshold": 2.0, "source": "PanelApp"},
     "column": "odds_ratio",
     "expect": 5.48,
     "tolerance": 0.05}
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
CLAIMS = Path(__file__).with_name("numbers.json")
MANUSCRIPT = Path(__file__).with_name("manuscript.md")


def _background_share(df: pd.DataFrame) -> float:
    """Share of positives among every gene outside the shared layer.

    The odds ratios in the paper are the shared layer against all other genes,
    not against genes of breadth zero, so the denominator has to be rebuilt
    from the dose-response table rather than read off one of its rows.
    """
    core = df[df.breadth == 6].iloc[0]
    hits = (df.fraction * df.n_genes).sum() - core.fraction * core.n_genes
    return float(hits / (df.n_genes.sum() - core.n_genes))


def _background_n(df: pd.DataFrame) -> float:
    return float(df.n_genes.sum() - df[df.breadth == 6].n_genes.iloc[0])


AGGREGATES = {
    "background_share": _background_share,
    "background_n": _background_n,
}


def lookup(claim: dict) -> float:
    df = pd.read_csv(RESULTS / claim["file"])
    for key, value in claim.get("where", {}).items():
        df = df[df[key] == value]

    if "aggregate" in claim:
        how = claim["aggregate"]
        if how not in AGGREGATES:
            raise LookupError(f"{claim['id']}: unknown aggregate {how!r}")
        return AGGREGATES[how](df)

    if len(df) != 1:
        raise LookupError(
            f"{claim['id']}: selector {claim.get('where')} matched {len(df)} rows "
            f"in {claim['file']}, expected exactly one"
        )
    return float(df.iloc[0][claim["column"]])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--verbose", action="store_true", help="print every comparison")
    args = ap.parse_args()

    claims = json.loads(CLAIMS.read_text(encoding="utf-8"))
    text = MANUSCRIPT.read_text(encoding="utf-8")

    failures: list[str] = []
    for claim in claims:
        actual = lookup(claim)
        expected = float(claim["expect"])
        tol = float(claim.get("tolerance", 0.01))
        ok = abs(actual - expected) <= tol
        in_text = claim["text"] in text
        if args.verbose or not (ok and in_text):
            mark = "ok  " if ok and in_text else "FAIL"
            print(f"{mark} {claim['id']:<34} manuscript={expected:<12.6g} "
                  f"results={actual:<12.6g} phrase={'found' if in_text else 'MISSING'}")
        if not ok:
            failures.append(f"{claim['id']}: manuscript says {expected}, "
                            f"{claim['file']} says {actual:.6g}")
        if not in_text:
            failures.append(f"{claim['id']}: phrase {claim['text']!r} not found in manuscript.md")

    print(f"\n{len(claims) - len(failures)} of {len(claims)} checks passed.")
    if failures:
        print("\nMismatches:")
        for f in failures:
            print(f"  - {f}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
