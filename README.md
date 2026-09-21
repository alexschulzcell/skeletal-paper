# How widely a gene acts across the skeleton predicts how widely its mutations cause disease

Analysis code, derived data and manuscript sources for the paper of the same
name (under submission at the *American Journal of Human Genetics*, Report
format).

For every gene we count at how many of six anatomically named bone-mineral-
density sites it acts — head, arms, legs, pelvis, femoral neck, lumbar spine.
That number, 0 to 6, is **breadth**. The paper shows four things about it:

1. **Breadth predicts Mendelian skeletal disease.** Genes that act at all six
   sites carry the disease burden; genes that act at one site do not. Three
   independently curated truth sides agree.
2. **Breadth predicts how far the disease extends.** Among genes that cause a
   Mendelian skeletal disease, breadth predicts the number of body regions the
   phenotype touches, and this is not an artefact of how well a gene is
   studied.
3. **It does not cross the organ boundary.** Breadth measured across seven
   subcortical brain volumes says nothing about skeletal phenotype extent,
   while skeletal breadth does. This resolves the published phenome-wide null
   result: the correspondence between common-variant and Mendelian pleiotropy
   is anatomically local.
4. **It reaches the clinic.** The result replicates in a cohort sharing no
   participants with the discovery cohort, and the shared layer carries
   fracture risk.

One gene property, two clinical ends: severely mutated it gives generalized
skeletal dysplasia in childhood, finely tuned it gives low bone density and
fracture risk in old age. Site-restricted genes do neither.

**Where the numbers live.** This README deliberately quotes none of them. Every
number the manuscript states is declared in
[`06_manuscript/numbers.json`](06_manuscript/numbers.json), with the table and
column it comes from, and `06_manuscript/check_numbers.py` recomputes each one
and exits non-zero on any disagreement. Prose that repeats a value goes stale
the first time the analysis is re-run; the self-test does not.

---

## Quickstart

```bash
git clone https://github.com/alexschulzcell/skeletal-paper.git
cd skeletal-paper

conda env create -f environment.yml     # Python 3.12, MAGMA 1.10
conda activate skeletal-breadth

make                    # prints the targets; runs nothing
make check              # are the inputs there? seconds
make verify             # do the manuscript and the tables agree? seconds
make all                # the whole thing, about eight hours
```

`make` on its own prints help and does nothing, because two of the stages take
hours. Run them in order, or run `make all` and leave it.

| command | what it does | time |
|---|---|---|
| `make data` | download every GWAS and reference file, then verify | 3–6 h, ~40 GB |
| `make genes` | harmonise, run MAGMA, assemble the gene matrix | **2–6 h** |
| `make layers` | breadth, loci, threshold sweep | < 1 min |
| `make truth` | the four truth sides | ~10 min |
| `make analysis` | `layers` + `truth` + the seven analyses | ~15 min |
| `make figures` | redraw the figures from the panel tables in `results/` | ~2 min |
| `make manuscript` | supplement, then the three submission PDFs | ~1 min |
| `make verify` | check every manuscript number against `results/` | seconds |
| `make all` | all of the above, in order | ~8 h |

Both long stages are **restartable**. Re-running skips any trait whose output
already exists, so an interrupted download or a killed MAGMA run costs only the
traits that had not finished.

**Hardware.** 8 GB RAM, four cores, and 60 GB of free disk (40 GB downloads,
20 GB MAGMA intermediates). More cores help only the MAGMA stage; raise it with
`make genes JOBS=8`, and give it 16 GB of RAM if you do. Nothing needs a GPU or
a cluster. Put the downloads elsewhere with
`make all SKELBREADTH_DATA=/scratch/skelbreadth`.

**Is it hung?** No. `make genes` writes progress to `results/magma_run.log`;
tail it. MAGMA prints nothing for long stretches by design.

Every script also runs on its own, from anywhere, and every script explains
itself:

```bash
python 04_analysis/02_symmetry_law.py --help
```

---

## The pipeline, stage by stage

Stages run in the order their directories are numbered, and every script reads
only what an earlier-numbered script wrote. There are no hidden dependencies
and no step that has to be run twice.

```
00_setup/          download everything, then refuse to continue if it is wrong
  01_download_sumstats.sh      GWAS summary statistics + md5sums
  02_download_references.sh    LD panel, gene locations, HPO, ClinVar, PanelApp, GO
  03_check_inputs.py           existence, size, checksum, column names
  04_prepare_atac_peaks.py     GSE252289 peak workbook -> per-element BED
  config.py                    every path, threshold and trait, in one place
  common.py                    shared loaders; no analysis lives here

01_gene_scores/    from summary statistics to one Z per gene per trait
  01_harmonise_sumstats.py     chr:pos -> rsID, p-values, N per trait
  02_run_magma.sh              THE LONG STEP - 2-6 h, parallel over traits
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
  01_enrichment.py             claim 1, three sources, four robustness levels
  02_symmetry_law.py           claim 2, with the ascertainment control
  03_cross_organ_control.py    claim 3, the 2x2
  04_confounders.py            length, publications, constraint, matched nulls
  05_replication_gefos.py      claim 4a, independent cohort
  06_clinical_endpoints.py     claim 4b, fracture and the other endpoints
  07_pathways_and_axis.py      mechanism: Wnt/ossification, regulatory density
  07b_recheck_regulatory_density.py   Figure 4B from GSE252289, no pipeline

05_figures/        the four figures, the graphical abstract, and two modules
  style.py                     colours, sizes, export, and the layout guards
  anatomy.py                   the body outline used by Figure 1A and the GA

06_manuscript/     manuscript, supplement, cover letter, and their checks
  manuscript.md                the text; the only place it is written
  numbers.json                 every numerical claim, with its source table
  check_numbers.py             recomputes each claim; the self-test that matters
  build_supplement.py          writes supplement.md from results/
  render_pdf.py                the three submission PDFs, with line numbers

results/           every table the paper cites            see results/README.md
figures/           every figure the paper shows           see figures/README.md
docs/              DATA_SOURCES.md, KNOWN_ISSUES.md
```

### How a number gets from the data into the paper

```
summary statistics -> MAGMA -> gene_z_matrix.tsv -> gene_layers.tsv
   -> 04_analysis/*.py -> results/<analysis>.tsv          [the analysis tables]
                       ~~~ gap, see below ~~~
   -> results/fig*.csv                                    [the panel tables]
   -> 05_figures/*.py -> figures/*.pdf, *.tif
   -> 06_manuscript/manuscript.md, checked against results/fig*.csv
```

The panel tables (`results/fig*.csv`) are the interface between the analysis
and everything downstream: the figures and the manuscript self-test read only
those, so a figure can be redrawn, and every number re-verified, without
re-running MAGMA. That half of the chain is complete and runs from a clean
checkout.

**The gap.** No script in this repository currently derives the panel tables
from the analysis tables. The `fig*.csv` files that ship here are the frozen
outputs of the analysis run the paper reports, and `make analysis` writes its
results beside them under different names rather than regenerating them. The
two sets also differ in content: the panel tables carry confidence intervals
and a logistic trend test that the analysis tables do not, so the missing step
is a real computation and not a rename. Until it is written,
**`make figures` redraws the published figures, but `make analysis` does not
reproduce their inputs.** This is issue 11 in
[`docs/KNOWN_ISSUES.md`](docs/KNOWN_ISSUES.md), which specifies what the
missing step has to do.

Before every submission, and after any change to the analysis, run

```bash
make verify
```

It is the self-test that matters: if it exits 0, the manuscript is consistent
with the tables in `results/`.

---

## Design rules this repository follows

- **Only what is in the paper.** Every script here produces a number or a
  figure the manuscript uses.
- **One definition, one place.** Breadth, the two layer thresholds, the six
  sites, the six body regions and the sweep grids are defined once, in
  `00_setup/config.py`. No script re-implements them.
- **Every number from a file.** Nothing is computed in a notebook. There are no
  notebooks.
- **No number typed into prose twice.** Documentation names files, not values.
  The manuscript is the one exception, and it is machine-checked.
- **No raw data in the repository.** `00_setup/` fetches it with checksums;
  `.gitignore` keeps it out.
- **No user paths.** Everything resolves from the repository root, overridable
  with `SKELBREADTH_ROOT`, `SKELBREADTH_DATA` and `MAGMA_BIN`.
- **A missing input fails at the top.** Every script checks its inputs before
  computing anything and names the file, the stage that produces it, and the
  command that would fix it.
- **Figures fail rather than ship broken.** `05_figures/style.py` refuses to
  write a figure whose canvas is not exactly the requested column width, whose
  ink falls outside that canvas, or that contains two pieces of overlapping
  text.

---

## What can and cannot be reproduced from a clean checkout

| stage | reproducible from public data? |
|---|---|
| the six DXA sites, heel BMD, fracture | yes, GWAS Catalog |
| GEFOS replication | yes, gefos.org |
| subcortical volumes, cross-organ control | yes, Oxford BIG40 |
| HPO, PanelApp, ClinVar, gnomAD, GO truth sides | yes |
| limb ATAC, regulatory density | yes, GEO GSE252289 — but **Figure 4B must be regenerated**, see issue 12 |
| **otosclerosis** | **no** — in-house harmonisation, not archived, and not reported in the paper. `06_clinical_endpoints.py` skips it with a stated message. |

---

## Limitations

We name them ourselves, in [`docs/KNOWN_ISSUES.md`](docs/KNOWN_ISSUES.md).
The ones a reader should know before reading a p-value:

- MAGMA reports no gene p below **5e-10**, so affected analyses are rank-based.
- The GEFOS sample sizes are **estimated** from standard errors, not read from
  a column.
- Redefining the layers inside GEFOS recovers only **19 %** of the same genes.
  The property replicates; the list does not, and the paper therefore reports
  breadth as a continuous property and never as a gene list.
- At **locus** level the direct layer contrast on fracture does not reach
  significance. The enrichment against background survives; the contrast does
  not.
- The MAGMA **annotation window** (35/10 kb) is a choice.
- Within the shared layer, fracture signal and Mendelian disease are **not**
  associated. The claim is about a layer, not about individual genes.
- The work was **exploratory, without pre-registration**.

---

## Licence and citation

Code is MIT ([`LICENSE`](LICENSE)). Derived data — everything in `results/`,
`figures/` and `docs/` — is CC BY 4.0 ([`LICENSE-DATA`](LICENSE-DATA)). The
upstream data are not redistributed here and carry their own terms; each is
named in [`docs/DATA_SOURCES.md`](docs/DATA_SOURCES.md).

To cite, see [`CITATION.cff`](CITATION.cff).
