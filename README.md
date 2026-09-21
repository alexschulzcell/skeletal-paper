# How widely a gene acts across the skeleton predicts how widely its mutations cause disease

Analysis code and derived result tables for the paper of the same name.

For every gene we count at how many of six anatomically named bone-mineral-density
sites it acts — head, arms, legs, pelvis, femoral neck, lumbar spine. That number,
0 to 6, is **breadth**. The paper shows four things about it:

1. **Breadth predicts Mendelian skeletal disease.** Genes acting at all six sites
   carry the disease burden; genes acting at one site do not. Three independently
   curated truth sides agree.
2. **Breadth predicts how far the disease extends**, and not because well-studied
   genes are better annotated.
3. **It does not cross the organ boundary.** Breadth across seven subcortical brain
   volumes says nothing about skeletal phenotype extent, while skeletal breadth does.
4. **It reaches the clinic.** The result replicates in a cohort sharing no
   participants with the discovery cohort, and the shared layer carries fracture risk.

This README quotes no numbers. Every value the paper reports is in a table under
[`results/`](results/), named in the docstring of the script that wrote it.

---

## Which script produced what

| Figure | Panel tables in `results/` | Written by |
|---|---|---|
| 1B–1D | `fig1b_*.csv`, `fig1c_*.csv`, `fig1d_*.csv` | `04_analysis/01_enrichment.py` |
| 2A | `fig2a_*.csv` | `04_analysis/02_symmetry_law.py` |
| 2B, 2C | `fig2b_*.csv`, `fig2c_*.csv` | `04_analysis/03_cross_organ_control.py` |
| 3A | `fig3a_*.csv` | `04_analysis/05_replication_gefos.py` |
| 3B, 3C | `fig3b_*.csv`, `fig3c_*.csv` | `04_analysis/06_clinical_endpoints.py` |
| 4A–4C | `fig4a_*.csv`, `fig4b_*.csv`, `fig4c_*.csv` | `04_analysis/07_pathways_and_axis.py` |

The confounder controls reported in the text come from `04_analysis/04_confounders.py`.
The figures themselves are drawn by `05_figures/fig1…fig4.py`, which read the panel
tables and compute nothing.

---

## The pipeline, stage by stage

Stages run in the order their directories are numbered, and every script reads only
what an earlier-numbered script wrote. Every script also runs on its own and explains
itself: `python 04_analysis/02_symmetry_law.py --help`.

```
00_setup/          fetch the inputs, then refuse to continue if they are wrong
  01_download_sumstats.sh      GWAS summary statistics + md5sums
  02_download_references.sh    LD panel, gene locations, HPO, ClinVar, PanelApp, GO
  03_check_inputs.py           existence, size, checksum, column names
  04_prepare_atac_peaks.py     GSE252289 peak workbook -> per-element BED
  config.py                    every path, threshold and trait, in one place
  common.py                    shared loaders; no analysis lives here

01_gene_scores/    from summary statistics to one Z per gene per trait
  01_harmonise_sumstats.py     chr:pos -> rsID, p-values, N per trait
  02_run_magma.sh              the long step - 2-6 h, parallel over traits
  03_assemble_gene_matrix.py   -> results/gene_z_matrix.tsv

02_layers/         the one number the paper is about
  01_define_breadth.py         breadth, peak, shared layer, site-specific layer
  02_locus_clumping.py         loci instead of genes, five merge distances
  03_threshold_sweep.py        the same layers at every threshold, and none

03_truth_sides/    four independent answers to "is this a disease gene?"
  01_build_hpo_regions.py      six body regions from the ontology
  02_build_panelapp.py         expert diagnostic panel, green genes
  03_build_clinvar.py          variant-level, with its own benign control
  04_build_brain_truth.py      for the cross-organ control only

04_analysis/       the seven analyses the paper reports
  01_enrichment.py             breadth against Mendelian disease
  02_symmetry_law.py           phenotype extent, with the ascertainment control
  03_cross_organ_control.py    the skeleton against the brain
  04_confounders.py            length, publications, constraint, matched nulls
  05_replication_gefos.py      independent cohort
  06_clinical_endpoints.py     fracture and the other endpoints
  07_pathways_and_axis.py      Wnt/ossification, regulatory density

05_figures/        the four figures and the graphical abstract
  style.py                     colours, sizes, export
  anatomy.py                   the body outline used by Figure 1A

results/           every table the paper cites            see results/README.md
figures/           every figure the paper shows           see figures/README.md
docs/              DATA_SOURCES.md, KNOWN_ISSUES.md
```

---

## Running it

```bash
conda env create -f environment.yml     # Python 3.12, MAGMA 1.10
conda activate skeletal-breadth

make                    # prints the targets; runs nothing
make check              # are the inputs there?        seconds
make all                # the whole thing              ~8 h
```

| command | what it does | time |
|---|---|---|
| `make data` | download every GWAS and reference file, then verify | 3–6 h, ~40 GB |
| `make genes` | harmonise, run MAGMA, assemble the gene matrix | **2–6 h** |
| `make layers` | breadth, loci, threshold sweep | < 1 min |
| `make truth` | the four truth sides | ~10 min |
| `make analysis` | `layers` + `truth` + the seven analyses | ~15 min |
| `make figures` | redraw the figures from the panel tables | ~2 min |

Both long stages are **restartable**: re-running skips any trait whose output already
exists. `make genes` writes progress to `results/magma_run.log`; MAGMA prints nothing
for long stretches by design.

**Hardware.** 8 GB RAM, four cores, 60 GB free disk. No GPU, no cluster. Raise
parallelism with `make genes JOBS=8`; put the downloads elsewhere with
`make all SKELBREADTH_DATA=/scratch/skelbreadth`.

---

## What can be reproduced from a clean checkout

| stage | public? |
|---|---|
| the six DXA sites, heel BMD, fracture | yes, GWAS Catalog |
| GEFOS replication | yes, gefos.org |
| subcortical volumes, cross-organ control | yes, Oxford BIG40 |
| HPO, PanelApp, ClinVar, gnomAD, GO truth sides | yes |
| limb ATAC, regulatory density | yes, GEO GSE252289 |

Every source is named, with version and licence, in
[`docs/DATA_SOURCES.md`](docs/DATA_SOURCES.md). The panel tables under `results/` are
the frozen outputs of the reported run; see [`results/README.md`](results/README.md).

---

## Limitations

Stated in full in [`docs/KNOWN_ISSUES.md`](docs/KNOWN_ISSUES.md). The ones to know
before reading a p-value:

- MAGMA reports no gene p below **5e-10**, so affected analyses are rank-based.
- The GEFOS sample sizes are **estimated** from standard errors, not read from a column.
- Redefining the layers inside GEFOS recovers only **19 %** of the same genes. The
  property replicates; the list does not, which is why breadth is reported as a
  continuous property and never as a gene list.
- At **locus** level the direct layer contrast on fracture does not reach significance.
- The MAGMA **annotation window** (35/10 kb) is a choice.
- Within the shared layer, fracture signal and Mendelian disease are **not** associated.
  The claim is about a layer, not about individual genes.
- The work was **exploratory, without pre-registration**.

---

## Licence and citation

Code is MIT ([`LICENSE`](LICENSE)). Derived data — everything in `results/`, `figures/`
and `docs/` — is CC BY 4.0 ([`LICENSE-DATA`](LICENSE-DATA)). Upstream data are not
redistributed here and carry their own terms. To cite, see [`CITATION.cff`](CITATION.cff).
