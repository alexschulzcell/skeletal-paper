"""Shared loaders for the skeletal-breadth companion repository.

This module holds the handful of operations that more than one stage needs, so
that they exist in exactly one place: reading MAGMA gene output, mapping Entrez
IDs to gene symbols, parsing the HPO ontology, and reading the two tables that
the earlier stages write into `results/`.

Nothing here performs an analysis. Every statistic of the paper is computed in
a numbered script under 02_layers/, 03_truth_sides/ or 04_analysis/.
"""

from __future__ import annotations

import collections
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as cfg  # noqa: E402


# --------------------------------------------------------------------------
# MAGMA gene-level output
# --------------------------------------------------------------------------

def gene_location_table() -> pd.DataFrame:
    """NCBI37.3 gene locations, indexed by Entrez ID.

    Columns: chr, start, end, strand, symbol.
    """
    path = cfg.require(cfg.REFERENCE / "NCBI37.3.gene.loc",
                       "MAGMA gene location file (NCBI37.3.gene.loc)")
    loc = pd.read_csv(path, sep="\t", header=None,
                      names=["id", "chr", "start", "end", "strand", "symbol"])
    return loc.drop_duplicates("id").set_index("id")


def read_magma_z(trait: str, column: str = "ZSTAT") -> pd.Series:
    """One column of MAGMA's `.genes.out` for `trait`, indexed by Entrez ID."""
    path = cfg.require(cfg.magma_out(trait), f"MAGMA gene results for trait '{trait}'")
    return pd.read_csv(path, sep=r"\s+", index_col=0)[column]


def magma_matrix(traits, extra=None) -> pd.DataFrame:
    """Gene x trait matrix of MAGMA ZSTAT, indexed by gene symbol.

    `traits` are required and rows missing any of them are dropped. `extra`
    traits are joined on top and may contain missing values, so that an
    endpoint that was not computed does not silently shrink the analysis set.
    Duplicate symbols are resolved by keeping the first Entrez ID, as MAGMA
    emits them, which is the convention used throughout the paper.
    """
    traits = list(traits)
    Z = pd.DataFrame({t: read_magma_z(t) for t in traits}).dropna()
    for t in (extra or []):
        Z[t] = read_magma_z(t)
    loc = gene_location_table().reindex(Z.index)
    Z["symbol"] = loc["symbol"].values
    Z["chr"] = loc["chr"].astype(str).values
    Z["gene_start"] = loc["start"].values
    Z["gene_end"] = loc["end"].values
    Z = Z.dropna(subset=["symbol", "chr"]).drop_duplicates("symbol").set_index("symbol")
    Z.index.name = "gene"
    return Z


# --------------------------------------------------------------------------
# Layer definitions - the one place where "core" and "site-specific" are made
# --------------------------------------------------------------------------

def breadth(Z: pd.DataFrame, sites=None, z=None) -> pd.Series:
    """Number of sites at which the gene passes the threshold (0-6)."""
    sites = list(sites or cfg.BMD_SITES)
    return (Z[sites] > (cfg.CORE_Z if z is None else z)).sum(axis=1)


def peak(Z: pd.DataFrame, sites=None) -> pd.Series:
    """Strongest single-site signal."""
    return Z[list(sites or cfg.BMD_SITES)].max(axis=1)


def core_genes(Z: pd.DataFrame, sites=None, z=None) -> set:
    """Shared core: above threshold at *every* site."""
    sites = list(sites or cfg.BMD_SITES)
    return set(Z[(Z[sites] > (cfg.CORE_Z if z is None else z)).all(axis=1)].index)


def specific_genes(Z: pd.DataFrame, sites=None, z0=None, delta=None,
                   exclude_core=True) -> set:
    """Site-specific layer: strong at one site, and clear of all the others."""
    sites = list(sites or cfg.BMD_SITES)
    z0 = cfg.SPECIFIC_Z if z0 is None else z0
    delta = cfg.SPECIFIC_DELTA if delta is None else delta
    out: set = set()
    for t in sites:
        others = [x for x in sites if x != t]
        out |= set(Z[(Z[t] > z0) & ((Z[t] - Z[others].max(axis=1)) > delta)].index)
    if exclude_core:
        out -= core_genes(Z, sites=sites)
    return out


def assign_layers(Z: pd.DataFrame, sites=None) -> pd.Series:
    """'core' / 'specific' / 'rest' for every gene in `Z`."""
    core = core_genes(Z, sites=sites)
    spec = specific_genes(Z, sites=sites)
    return pd.Series(
        np.where(Z.index.isin(core), "core",
                 np.where(Z.index.isin(spec), "specific", "rest")),
        index=Z.index, name="layer")


# --------------------------------------------------------------------------
# Tables written by earlier stages
# --------------------------------------------------------------------------

def load_gene_matrix() -> pd.DataFrame:
    """results/gene_z_matrix.tsv, from 01_gene_scores/03_assemble_gene_matrix.py."""
    path = cfg.require(cfg.RESULTS / "gene_z_matrix.tsv",
                       "gene x trait Z matrix")
    return pd.read_csv(path, sep="\t", index_col=0)


def load_layers() -> pd.DataFrame:
    """results/gene_layers.tsv, from 02_layers/01_define_breadth.py."""
    path = cfg.require(cfg.RESULTS / "gene_layers.tsv",
                       "gene layer assignment (breadth, peak, layer)")
    return pd.read_csv(path, sep="\t", index_col=0)


def load_truth(name: str) -> pd.DataFrame:
    """One of the truth-side tables written by 03_truth_sides/."""
    path = cfg.require(cfg.RESULTS / f"truth_{name}.tsv",
                       f"truth-side table '{name}'")
    return pd.read_csv(path, sep="\t", index_col=0)


# --------------------------------------------------------------------------
# HPO ontology
# --------------------------------------------------------------------------

def hpo_descendants(roots) -> dict:
    """Map every HPO root term to its full descendant set (root included).

    Reads `hp.obo`. `roots` is a mapping name -> list of HP:xxxxxxx terms, as in
    config.HPO_SKELETAL_REGIONS; the return value has the same keys.
    """
    path = cfg.require(cfg.REFERENCE / "hp.obo", "Human Phenotype Ontology (hp.obo)")
    parents = collections.defaultdict(list)
    current = None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line == "[Term]":
                current = None
            elif line.startswith("id: HP:"):
                current = line[4:]
            elif line.startswith("is_a:") and current:
                parents[current].append(line[6:16])
    children = collections.defaultdict(list)
    for child, ps in parents.items():
        for p in ps:
            children[p].append(child)

    def descend(root):
        seen, stack = set(), [root]
        while stack:
            node = stack.pop()
            if node in seen:
                continue
            seen.add(node)
            stack.extend(children.get(node, []))
        return seen

    return {name: set().union(*[descend(r) for r in terms])
            for name, terms in roots.items()}


def gene_to_hpo_terms() -> dict:
    """Gene symbol -> set of HPO terms, via genes_to_disease and phenotype.hpoa.

    Only phenotype ("P") rows of phenotype.hpoa are used, so that the mode of
    inheritance and clinical-course annotations do not inflate the term count
    that the ascertainment control conditions on.
    """
    hpoa = cfg.require(cfg.REFERENCE / "phenotype.hpoa", "HPO disease annotations")
    g2d_path = cfg.require(cfg.REFERENCE / "genes_to_disease.txt",
                           "HPO gene-to-disease mapping")

    disease_terms = collections.defaultdict(set)
    with open(hpoa, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(("#", "database_id")):
                continue
            f = line.split("\t")
            if len(f) > 10 and f[10] == "P":
                disease_terms[f[0]].add(f[3])

    gene_diseases = collections.defaultdict(set)
    with open(g2d_path, encoding="utf-8") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > 3:
                gene_diseases[f[1]].add(f[3])

    return {g: set().union(*[disease_terms.get(d, set()) for d in ds]) if ds else set()
            for g, ds in gene_diseases.items()}


# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------

def write_result(df: pd.DataFrame, name: str, index: bool = True) -> Path:
    """Write one table into results/ and say so on stdout."""
    cfg.ensure_results_dir()
    path = cfg.RESULTS / name
    df.to_csv(path, sep="\t", index=index)
    print(f"  -> {path.relative_to(cfg.ROOT)}  ({len(df)} rows)")
    return path


def banner(title: str) -> None:
    print("\n" + "=" * 74)
    print(title)
    print("=" * 74)
