# `results/`

Every number the manuscript states comes from a file in this directory, and
from nowhere else. No number is computed in a notebook, typed in by hand, or
read off a figure. `06_manuscript/check_numbers.py` enforces exactly that: it
recomputes each claim from the table it came from and exits non-zero on the
first mismatch.

Two families of file live here, and only one of them is currently regenerated
by the pipeline. See issue 11 in `docs/KNOWN_ISSUES.md`.

## 1 · Analysis tables — written by the pipeline

Each is a tab-separated table written by exactly one script. To find out which,
search the repository for the filename: it appears in the docstring of the
script that writes it, under "Writes". The top-level `README.md` maps every
number in the paper to its script and its table.

| written by | files |
|---|---|
| `00_setup/` | `input_check.txt` |
| `01_gene_scores/` | `harmonisation_report.tsv`, `gene_z_matrix.tsv`, `gene_coverage.tsv`, `magma_run.log` |
| `02_layers/` | `gene_layers.tsv`, `layer_summary.tsv`, `breadth_distribution.tsv`, `loci.tsv`, `locus_table_*.tsv`, `locus_summary.tsv`, `locus_endpoint_contrast.tsv`, `threshold_layers.tsv`, `threshold_layer_sizes.tsv`, `continuous_axes.tsv` |
| `03_truth_sides/` | `truth_hpo.tsv`, `truth_panelapp.tsv`, `truth_clinvar.tsv`, `truth_brain.tsv`, `hpo_region_sizes.tsv`, `clinvar_provenance.txt` |
| `04_analysis/` | `enrichment_*.tsv`, `symmetry_*.tsv`, `cross_organ_*.tsv`, `confounder_*.tsv`, `replication_*.tsv`, `endpoints_*.tsv`, `pathway_*.tsv`, `regulatory_density.tsv`, `maturation_axis*.tsv` |

Only tables the paper cites belong here. Intermediate products that no figure
and no sentence depends on are not written out.

## 2 · Panel tables — one file per figure panel

The `fig*_*.csv` files are the figure layer's inputs: one comma-separated table
per panel, with English column names, holding exactly the values that panel
draws. The scripts in `05_figures/` read only these and compute nothing, so a
figure can be redrawn without re-running the analysis, and a number in a figure
can be checked by opening one small file.

| panel table | backed by the analysis table |
|---|---|
| `fig1b_dose_response.csv`, `fig1b_trend_tests.csv` | `enrichment_dose_response.tsv` |
| `fig1c_threshold_sweep.csv` | `enrichment_threshold_sweep.tsv` |
| `fig1d_locus_sweep.csv` | `enrichment_locus_sweep.tsv` |
| `fig2a_symmetry_law.csv`, `fig2a_symmetry_stats.csv`, `fig2a_symmetry_regression.csv` | `symmetry_by_breadth.tsv`, `symmetry_models.tsv` |
| `fig2b_cross_organ.csv` | `cross_organ_matrix.tsv` |
| `fig2c_trait_coherence.csv` | `cross_organ_correlation.tsv` |
| `fig3a_gefos_replication.csv`, `fig3a_gefos_stats.csv`, `fig3a_gefos_conditional.csv` | `replication_dose_response.tsv`, `replication_layers.tsv`, `replication_models.tsv` |
| `fig3b_fracture.csv` | `endpoints_dose_response.tsv`, `endpoints_by_layer.tsv` |
| `fig3c_layer_contrast.csv`, `fig3c_conditional.csv` | `endpoints_models.tsv` |
| `fig4a_pathways.csv` | `pathway_enrichment.tsv` |
| `fig4b_regulatory_density.csv` | `regulatory_density.tsv` |
| `fig4c_shared_layer_genes.csv` | `gene_layers.tsv`, `pathway_core_members.tsv` |

If a panel table and the analysis table behind it disagree, the analysis table
is right and the panel table is stale.

**The link is not automated.** No script derives the `fig*.csv` files from the
analysis tables; the ones shipped here are the frozen outputs of the run the
paper reports. The panel tables also carry quantities the analysis tables do
not — Jeffreys and Wald intervals, a logistic trend test, and a layer contrast
that excludes the other layer from the comparison group — so the missing step
is a computation, not a rename. `docs/KNOWN_ISSUES.md` issue 11 specifies it.

---

The tables are licensed CC BY 4.0; see `LICENSE-DATA`.
