#!/usr/bin/env python3
"""Fill the submission folder, and put nothing else in it.

    python 06_manuscript/assemble_submission.py
    python 06_manuscript/assemble_submission.py --dest /path/to/submission

What lands there is exactly what the journal is sent: the manuscript, the
supplemental information and the cover letter as PDFs, and the five figures as
TIFFs. Working files, checklists, notes and the vector copies of the figures
stay in this repository, where they belong.

The folder is emptied first. Anything a person dropped in by hand is reported
and removed, because a stray file in a submission folder is how the wrong
version reaches an editor.

File names carry no spaces. Submission systems mangle them, and Cell Press
asks for figures named Figure1, Figure2 and so on.
"""
from __future__ import annotations

import argparse
import hashlib
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PDF = ROOT / "06_manuscript" / "pdf"
FIG = ROOT / "figures"

# destination name -> source path
MANIFEST: dict[str, Path] = {
    "Manuscript.pdf": PDF / "Manuscript.pdf",
    "SupplementalInformation.pdf": PDF / "Supplemental Information.pdf",
    "CoverLetter.pdf": PDF / "Cover Letter.pdf",
    "Figure1.tif": FIG / "figure1.tif",
    "Figure2.tif": FIG / "figure2.tif",
    "Figure3.tif": FIG / "figure3.tif",
    "Figure4.tif": FIG / "figure4.tif",
    "GraphicalAbstract.tif": FIG / "graphical_abstract.tif",
}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:12]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dest", type=Path, default=ROOT.parent / "submission",
                    help="submission folder (default: ../submission)")
    args = ap.parse_args()
    dest: Path = args.dest

    missing = [n for n, src in MANIFEST.items() if not src.exists()]
    if missing:
        print("Nothing copied. These build products do not exist yet:")
        for n in missing:
            print(f"  {n:32s} <- {MANIFEST[n].relative_to(ROOT)}")
        print("\nRun:  make figures && make manuscript")
        return 1

    dest.mkdir(parents=True, exist_ok=True)
    removed = []
    for old in sorted(dest.iterdir()):
        if old.is_dir():
            shutil.rmtree(old)
            removed.append(old.name + "/")
        else:
            old.unlink()
            removed.append(old.name)
    if removed:
        print(f"emptied {dest}:")
        for r in removed:
            print(f"  removed {r}")

    print(f"\nwriting {dest}:")
    total = 0
    for name, src in MANIFEST.items():
        shutil.copy2(src, dest / name)
        size = (dest / name).stat().st_size
        total += size
        print(f"  {name:32s} {size / 1e6:7.2f} MB  sha256:{digest(dest / name)}")
    print(f"  {'':32s} {total / 1e6:7.2f} MB total, {len(MANIFEST)} files")

    over = [n for n in MANIFEST if (dest / n).stat().st_size > 15e6]
    if over:
        print("\nWARNING: over the 15 MB per-file limit: " + ", ".join(over))
    return 0


if __name__ == "__main__":
    sys.exit(main())
