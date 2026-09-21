#!/usr/bin/env python3
"""Verify that every input the pipeline needs is present, readable and intact.

What it does
    Walks the list of files the two download scripts produce and reports, for
    each, whether it exists, whether it is non-empty, whether its md5 matches
    the recorded checksum, and - for the summary statistics - whether the
    columns the harmonisation step needs are actually there. It then prints a
    verdict and exits non-zero if anything required is missing, so that
    `make all` stops here rather than four hours into MAGMA.

Reads
    $SKELBREADTH_DATA/sumstats/*.tsv.gz, and md5sums.txt beside them
    $SKELBREADTH_DATA/reference/*

Writes
    results/input_check.txt   the same report, for the record

Usage
    python 00_setup/03_check_inputs.py [--quick] [--no-checksums]
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as cfg  # noqa: E402

# (filename, required?, human description)
SUMSTATS = [
    ("head.tsv.gz",   True,  "DXA BMD, head (GCST90568448)"),
    ("arms.tsv.gz",   True,  "DXA BMD, arms (GCST90568449)"),
    ("legs.tsv.gz",   True,  "DXA BMD, legs (GCST90568450)"),
    ("pelvis.tsv.gz", True,  "DXA BMD, pelvis (GCST90568451)"),
    ("fn.tsv.gz",     True,  "DXA BMD, femoral neck (GCST90568452)"),
    ("ls.tsv.gz",     True,  "DXA BMD, lumbar spine (GCST90568453)"),
    ("heel.tsv.gz",   True,  "Heel BMD, ultrasound (GCST90568454)"),
    ("frak.tsv.gz",   True,  "Fracture (GCST006980)"),
    ("height.tsv.gz", True,  "Standing height (Yengo 2022)"),
    ("prop.tsv.gz",   True,  "Sitting-height ratio (GCST90728588)"),
    ("fn2stu.tsv.gz", True,  "GEFOS femoral neck (Zheng 2015)"),
    ("ls2stu.tsv.gz", True,  "GEFOS lumbar spine (Zheng 2015)"),
    ("fa2stu.tsv.gz", True,  "GEFOS forearm (Zheng 2015)"),
] + [(f"{r}.tsv.gz", True, f"Subcortical volume, {r} (Oxford BIG40)")
     for r in cfg.BRAIN_REGIONS]

REFERENCE = [
    ("g1000_eur.bim",           True,  "1000 Genomes EUR LD reference (MAGMA)"),
    ("g1000_eur.bed",           True,  "1000 Genomes EUR LD reference (MAGMA)"),
    ("g1000_eur.fam",           True,  "1000 Genomes EUR LD reference (MAGMA)"),
    ("NCBI37.3.gene.loc",       True,  "Gene locations, build 37"),
    ("hp.obo",                  True,  "Human Phenotype Ontology"),
    ("phenotype.hpoa",          True,  "HPO disease annotations"),
    ("genes_to_disease.txt",    True,  "HPO gene-to-disease mapping"),
    ("clinvar.vcf.gz",          True,  "ClinVar variant records (GRCh37)"),
    ("panelapp_pa309.json",     True,  "PanelApp PA309, skeletal dysplasia"),
    ("gnomad_constraint.txt.bgz", True, "gnomAD v2.1.1 LOEUF"),
    ("goa_human.gaf.gz",        False, "GO annotations - needed by 04_analysis/07"),
    ("go-basic.obo",            False, "GO ontology - needed by 04_analysis/07"),
    ("gene2pubmed.gz",          False, "Publication counts - needed by 04_analysis/04"),
]

# Columns the harmonisation step will look for, by family of source.
EXPECTED_COLUMNS = {
    "ebi":   {"p_value", "chromosome", "base_pair_location"},
    "gefos": {"p-value", "chromosome", "position"},
    "big40": {"pval"},
}


def md5(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def header_of(path: Path) -> list[str]:
    """First line of a gzipped table, split on whitespace."""
    try:
        with gzip.open(path, "rt", errors="ignore") as fh:
            return fh.readline().strip().replace(",", "\t").split()
    except OSError:
        return []


def check_group(name, directory, entries, do_checksums, quick, report):
    def say(line):
        print(line)
        report.append(line)

    say(f"\n=== {name} ===")
    say(f"    {directory}")
    recorded = {}
    sums = directory / "md5sums.txt"
    if do_checksums and sums.exists():
        for line in sums.read_text().splitlines():
            parts = line.split()
            if len(parts) == 2:
                recorded[Path(parts[1]).name] = parts[0]

    missing_required, warnings = [], []
    for filename, required, description in entries:
        path = directory / filename
        if not path.exists():
            mark, note = ("MISSING" if required else "absent  "), description
            (missing_required if required else warnings).append((filename, description))
            say(f"  [{mark:7s}] {filename:24s} {note}")
            continue
        size = path.stat().st_size
        if size == 0:
            missing_required.append((filename, "file is empty"))
            say(f"  [EMPTY  ] {filename:24s} {description}")
            continue
        note = f"{size / 1e6:9.1f} MB"
        if do_checksums and filename in recorded and not quick:
            got = md5(path)
            if got != recorded[filename]:
                missing_required.append((filename, "md5 mismatch"))
                say(f"  [BAD MD5] {filename:24s} {note}  expected {recorded[filename]}")
                continue
            note += "  md5 ok"
        if filename.endswith(".tsv.gz") and not quick:
            cols = set(header_of(path))
            if cols and not any(exp & cols for exp in EXPECTED_COLUMNS.values()):
                warnings.append((filename, f"unrecognised columns: {sorted(cols)[:6]}"))
                note += "  (unrecognised header)"
        say(f"  [ok     ] {filename:24s} {note}  {description}")
    return missing_required, warnings


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--quick", action="store_true",
                    help="existence and size only; skip md5 and header checks")
    ap.add_argument("--no-checksums", action="store_true",
                    help="skip md5 verification but still read headers")
    args = ap.parse_args()

    report: list[str] = ["Input check for the skeletal-breadth repository"]
    missing, warnings = [], []

    for name, directory, entries in [
        ("Summary statistics", cfg.SUMSTATS, SUMSTATS),
        ("Reference files", cfg.REFERENCE, REFERENCE),
    ]:
        if not directory.exists():
            line = f"\n=== {name} ===\n  DIRECTORY DOES NOT EXIST: {directory}"
            print(line)
            report.append(line)
            missing.append((str(directory), "directory missing"))
            continue
        m, w = check_group(name, directory, entries,
                           not args.no_checksums, args.quick, report)
        missing += m
        warnings += w

    def say(line):
        print(line)
        report.append(line)

    say("\n" + "=" * 74)
    if warnings:
        say(f"{len(warnings)} optional input(s) absent or unusual:")
        for f, why in warnings:
            say(f"  - {f}: {why}")
        say("  The stages that need them will say so and skip; see the README table.")
    if missing:
        say(f"{len(missing)} REQUIRED input(s) missing or corrupt:")
        for f, why in missing:
            say(f"  - {f}: {why}")
        say("\n  Run  bash 00_setup/01_download_sumstats.sh")
        say("  and  bash 00_setup/02_download_references.sh")
        say("  or see docs/DATA_SOURCES.md to obtain them another way.")
    else:
        say("All required inputs present. The pipeline can run.")
    say("=" * 74)

    cfg.ensure_results_dir()
    (cfg.RESULTS / "input_check.txt").write_text("\n".join(report) + "\n",
                                                 encoding="utf-8")
    print(f"\n  -> results/input_check.txt")
    return 1 if missing else 0


if __name__ == "__main__":
    sys.exit(main())
