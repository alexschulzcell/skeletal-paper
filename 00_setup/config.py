"""Central configuration for the skeletal-breadth companion repository.

Every script in this repository imports this module instead of hard-coding a
path. Nothing here depends on the machine it runs on: the repository root is
found from this file's own location, and the two directories that hold data
downloaded from public archives can be redirected with environment variables.

Environment variables
---------------------
SKELBREADTH_ROOT    override the repository root (default: parent of 00_setup/)
SKELBREADTH_DATA    where 00_setup/ downloads GWAS summary statistics and
                    reference files to (default: <root>/data_raw)
MAGMA_BIN           path to the MAGMA executable (default: `magma` on PATH)

Layer definitions
-----------------
The two layers of the paper are defined here, once, and nowhere else:

    shared core     MAGMA Z > CORE_Z at **all six** DXA sites
    site-specific   MAGMA Z > SPECIFIC_Z at one site and at least
                    SPECIFIC_DELTA above the best of the other five

`breadth` is the count of sites with Z > CORE_Z (0-6); `peak` is the maximum Z
across the six sites. Both are reported as continuous properties in the paper
(see docs/KNOWN_ISSUES.md on why the gene *list* is not the unit of claim).
"""

from __future__ import annotations

import os
from pathlib import Path

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------

ROOT = Path(os.environ.get("SKELBREADTH_ROOT", Path(__file__).resolve().parent.parent))
DATA = Path(os.environ.get("SKELBREADTH_DATA", ROOT / "data_raw"))

SUMSTATS = DATA / "sumstats"        # raw GWAS summary statistics, as downloaded
REFERENCE = DATA / "reference"      # 1000G, gene locations, HPO, ClinVar, PanelApp
MAGMA_DIR = DATA / "magma"          # MAGMA inputs and per-trait .genes.out
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"

MAGMA_BIN = os.environ.get("MAGMA_BIN", "magma")

# --------------------------------------------------------------------------
# Traits
# --------------------------------------------------------------------------

#: The six anatomically named DXA bone-mineral-density sites. Order is fixed and
#: is the order used in every table and figure of the paper.
BMD_SITES = ["head", "arms", "legs", "pelvis", "fn", "ls"]

BMD_SITE_LABELS = {
    "head": "Head",
    "arms": "Arms",
    "legs": "Legs",
    "pelvis": "Pelvis",
    "fn": "Femoral neck",
    "ls": "Lumbar spine",
}

#: Clinical and anthropometric endpoints, kept separate from the six sites that
#: define breadth so that no endpoint can enter its own predictor.
ENDPOINTS = {
    "frak": "Fracture",
    "heel": "Heel BMD (ultrasound)",
    "height": "Standing height",
    "prop": "Sitting-height ratio",
}

#: Replication cohort (GEFOS/UK10K, Zheng 2015). No UK Biobank participants.
GEFOS_TRAITS = {
    "fn2stu": "Femoral neck",
    "ls2stu": "Lumbar spine",
    "fa2stu": "Forearm",
}

#: Seven subcortical brain volumes (Oxford BIG40). Used *only* for the
#: cross-organ control in 04_analysis/03_cross_organ_control.py.
BRAIN_REGIONS = [
    "thalamus", "caudate", "putamen", "pallidum",
    "hippocampus", "amygdala", "accumbens",
]

#: Sample sizes, for MAGMA's `N=` argument and for docs/DATA_SOURCES.md.
TRAIT_N = {
    "head": 31986, "arms": 31986, "legs": 31986,
    "pelvis": 31986, "fn": 31986, "ls": 31986,
    "heel": 426824,
    "frak": 426795,
    "height": 1232747,
    "prop": 473511,
    "fn2stu": 34806, "ls2stu": 25759, "fa2stu": 8112,
    **{r: 33211 for r in ["thalamus", "caudate", "putamen", "pallidum",
                          "hippocampus", "amygdala", "accumbens"]},
}

# --------------------------------------------------------------------------
# Layer definitions
# --------------------------------------------------------------------------

CORE_Z = 2.0            # Z threshold that defines "acts at this site"
SPECIFIC_Z = 4.0        # Z threshold for the site-specific layer
SPECIFIC_DELTA = 2.0    # margin over the best other site
THRESHOLD_SWEEP = [2.0, 2.5, 3.0, 3.5, 4.0]     # 02_layers/03_threshold_sweep.py
CLUMP_GAPS = [0, 100_000, 250_000, 500_000, 1_000_000]  # 02_layers/02_locus_clumping.py
CLUMP_PRIMARY_GAP = 250_000

MAGMA_WINDOW = "35,10"  # kb upstream, kb downstream; see docs/KNOWN_ISSUES.md
MAGMA_P_FLOOR = 5e-10   # MAGMA's gene-p floor; see docs/KNOWN_ISSUES.md

# --------------------------------------------------------------------------
# Truth sides
# --------------------------------------------------------------------------

#: Six anatomical body regions on the disease side, mirroring the six GWAS
#: sites. Each entry is an HPO root whose full descendant set is taken.
HPO_SKELETAL_REGIONS = {
    "Skull":        ["HP:0000929"],
    "Spine":        ["HP:0000925"],
    "Thorax":       ["HP:0000765"],
    "Pelvis":       ["HP:0002644", "HP:0008839"],
    "Long bones":   ["HP:0011314"],
    "Hands/feet":   ["HP:0001155", "HP:0001760", "HP:0011297", "HP:0001780"],
}

#: Brain counterpart, for the cross-organ control only.
HPO_BRAIN_REGIONS = {
    "Cortex":         ["HP:0002060"],
    "Basal ganglia":  ["HP:0002134"],
    "Cerebellum":     ["HP:0001317"],
    "Brainstem":      ["HP:0002363"],
    "Corpus callosum": ["HP:0001273"],
    "Ventricles":     ["HP:0002119"],
    "White matter":   ["HP:0002500"],
    "Hippocampus":    ["HP:0025100"],
}

PANELAPP_PANEL_ID = 309         # "Skeletal dysplasia", confidence level 3 only
PANELAPP_CONFIDENCE = "3"

#: Regular expression matching a skeletal diagnosis in a ClinVar CLNDN field.
CLINVAR_SKELETAL_PATTERN = (
    r"dysplas|skelet|bone|chondro|osteo|craniosynostos|achondro|dwarf|"
    r"spondyl|brachydactyl|polydactyl|syndactyl|exostos|osteogenesis|"
    r"rickets|hypophosphat|cleidocranial|metaphys|epiphys|diaphys"
)

# --------------------------------------------------------------------------
# Reproducibility
# --------------------------------------------------------------------------

SEED = 2026
N_PERMUTATIONS = 5000   # matched-null draws in 04_analysis/04_confounders.py
                        # and 04_analysis/07_pathways_and_axis.py


# --------------------------------------------------------------------------
# Small helpers used by every script
# --------------------------------------------------------------------------

def require(path, what: str):
    """Abort with a readable message if an input file is missing.

    Every script calls this before touching an input, so that a missing
    download fails at the top with an instruction rather than deep inside
    pandas with a traceback.
    """
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        raise SystemExit(
            f"\nMISSING INPUT: {what}\n"
            f"  expected at: {path}\n"
            f"  Run the setup stage first:   make data\n"
            f"  See docs/DATA_SOURCES.md for what this file is and where it comes from.\n"
        )
    return path


def magma_out(trait: str) -> Path:
    """Path of MAGMA's gene-level output for one trait."""
    return MAGMA_DIR / f"{trait}.genes.out"


def ensure_results_dir() -> Path:
    RESULTS.mkdir(parents=True, exist_ok=True)
    return RESULTS
