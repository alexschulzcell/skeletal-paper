#!/usr/bin/env python3
"""Build the second truth side: the Genomics England skeletal dysplasia panel.

What it does
    Reads the PanelApp record for panel 309, "Skeletal dysplasia", and keeps
    the genes rated at confidence level 3 - green, diagnostic-grade, expert
    reviewed. This is the strictest of the paper's three curations and the one
    with the cleanest provenance: a clinical panel maintained for diagnostic
    use, not a text-mined list.

    It is deliberately the smallest of the three. Using all confidence levels
    would roughly double it and weaken the contrast by mixing in amber and red
    entries, which is the opposite of what a positive control should do.

Reads
    $SKELBREADTH_DATA/reference/panelapp_pa309.json

Writes
    results/truth_panelapp.tsv    gene, confidence, mode of inheritance,
                                  entity type

Runtime
    Seconds.

Usage
    python 03_truth_sides/02_build_panelapp.py
    python 03_truth_sides/02_build_panelapp.py --confidence 2 3   # sensitivity
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "00_setup"))
import common  # noqa: E402
import config as cfg  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--confidence", nargs="*", default=[cfg.PANELAPP_CONFIDENCE],
                    help="confidence levels to keep (default: 3, green only)")
    ap.add_argument("--out", default="truth_panelapp.tsv")
    args = ap.parse_args()

    common.banner("Truth side 2 of 4: PanelApp PA309, skeletal dysplasia")

    path = cfg.require(cfg.REFERENCE / "panelapp_pa309.json",
                       "PanelApp panel 309 (skeletal dysplasia)")
    with open(path, encoding="utf-8") as fh:
        panel = json.load(fh)

    print(f"  panel      : {panel.get('name', 'PA309')}")
    print(f"  version    : {panel.get('version', 'unknown')}")
    print(f"  last update: {panel.get('version_created', 'unknown')}")

    rows = []
    for entry in panel.get("genes", []):
        symbol = (entry.get("gene_data") or {}).get("gene_symbol")
        if not symbol:
            continue
        rows.append({
            "gene": symbol,
            "confidence": str(entry.get("confidence_level", "")),
            "mode_of_inheritance": entry.get("mode_of_inheritance", ""),
            "entity_type": entry.get("entity_type", "gene"),
        })
    d = pd.DataFrame(rows).drop_duplicates("gene").set_index("gene")

    print(f"\n  entries on the panel at any confidence: {len(d):,}")
    counts = d.confidence.value_counts().sort_index()
    for level, n in counts.items():
        label = {"3": "green, diagnostic grade", "2": "amber", "1": "red"}.get(
            level, "unrated")
        print(f"    confidence {level or '-':2s}: {n:5,d}   {label}")

    keep = d[d.confidence.isin([str(c) for c in args.confidence])].copy()
    keep["on_panel"] = 1
    print(f"\n  kept at confidence {', '.join(map(str, args.confidence))}: "
          f"{len(keep):,} genes")

    if len(keep) < 100:
        print("\n  WARNING: fewer genes than expected for PA309. The API layout may")
        print("  have changed; check the JSON before using this as a truth side.")

    if "mode_of_inheritance" in keep:
        moi = keep.mode_of_inheritance.str.extract(
            r"(BIALLELIC|MONOALLELIC|X-LINKED|MITOCHONDRIAL)", expand=False)
        print("\n=== Mode of inheritance, for the record ===")
        for mode, n in moi.value_counts().items():
            print(f"  {mode:16s}{n:5,d}")

    common.write_result(keep, args.out)
    print("\nNow run:  python 03_truth_sides/03_build_clinvar.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
