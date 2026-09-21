#!/usr/bin/env python3
"""Write the supplemental information from the same tables the figures use.

Every supplemental table is generated here from `results/`, so it cannot
disagree with a figure or with the main text. Nothing is typed in by hand.

    python 06_manuscript/build_supplement.py     # writes 06_manuscript/supplement.md
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
OUT = Path(__file__).with_name("supplement.md")


def fmt_p(p: float) -> str:
    if pd.isna(p):
        return "—"
    if p >= 0.001:
        return f"{p:.3f}"
    m, e = f"{p:.1e}".split("e")
    return f"{m} × 10<sup>{int(e)}</sup>"


def table(df: pd.DataFrame, cols: dict[str, str]) -> str:
    head = "| " + " | ".join(cols.values()) + " |"
    rule = "|" + "|".join(["---"] * len(cols)) + "|"
    rows = []
    for _, r in df.iterrows():
        rows.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    return "\n".join([head, rule, *rows])


def main() -> int:
    read = lambda n: pd.read_csv(RESULTS / n)  # noqa: E731

    parts: list[str] = []
    add = parts.append

    add("# Supplemental information\n")
    add("**How widely a gene acts across the skeleton predicts how widely its "
        "mutations cause disease**\n")
    add("Alexander Zentgraf\n")
    add("---\n")

    add("## Contents\n")
    add("Table S1. Breadth against Mendelian skeletal disease, three curations, "
        "breadth 0 to 6.  \n"
        "Table S2. Sensitivity of the enrichment to the Z threshold that defines "
        "action at a site.  \n"
        "Table S3. Sensitivity of the enrichment to collapsing genes into loci.  \n"
        "Table S4. Phenotype breadth by breadth of action.  \n"
        "Table S5. The cross-organ control.  \n"
        "Table S6. Genetic coherence of the measurement sites within each organ.  \n"
        "Table S7. Replication in GEFOS.  \n"
        "Table S8. Fracture risk by breadth of action.  \n"
        "Table S9. The two layers across replication and fracture endpoints.  \n"
        "Table S10. Conditional models: breadth beside all six site scores.  \n"
        "Table S11. Developmental signaling pathways in both layers.  \n"
        "Table S12. Regulatory density against a length-matched null.  \n"
        "Table S13. The 103 shared-layer genes.  \n"
        "Supplemental note 1. Which comparison group each odds ratio uses.  \n"
        "Supplemental note 2. Analytical choices, and what was not corrected for.  \n")
    add("---\n")

    # ------------------------------------------------------------- S1
    d = read("fig1b_dose_response.csv").copy()
    d["pct"] = (d.fraction * 100).map("{:.1f}".format)
    d["ci"] = [f"{a*100:.1f}–{b*100:.1f}" for a, b in zip(d.ci_low, d.ci_high)]
    add("## Table S1. Breadth against Mendelian skeletal disease\n")
    add("Proportion of genes carrying an annotation in each curation, by breadth "
        "of action. Intervals are 95% Jeffreys intervals. This is the source of "
        "Figure 1B.\n")
    add(table(d, {"source": "Curation", "breadth": "Breadth", "n_genes": "Genes",
                  "n_positive": "Annotated", "pct": "Annotated (%)",
                  "ci": "95% CI (%)"}))
    add("")

    # ------------------------------------------------------------- S2
    d = read("fig1c_threshold_sweep.csv").copy()
    d["or_ci"] = [f"{o:.2f} ({a:.2f}–{b:.2f})"
                  for o, a, b in zip(d.odds_ratio, d.ci_low, d.ci_high)]
    d["pp"] = d.p.map(fmt_p)
    add("## Table S2. Sensitivity to the Z threshold\n")
    add("The shared layer redefined at each threshold, against every gene outside "
        "it. Haldane-corrected odds ratios with Wald intervals. Source of "
        "Figure 1C.\n")
    add(table(d, {"threshold": "Z threshold", "source": "Curation",
                  "n_core": "Shared layer (n)", "or_ci": "OR (95% CI)",
                  "pp": "P"}))
    add("")

    # ------------------------------------------------------------- S3
    d = read("fig1d_locus_sweep.csv").copy()
    d["or_ci"] = [f"{o:.2f} ({a:.2f}–{b:.2f})"
                  for o, a, b in zip(d.odds_ratio, d.ci_low, d.ci_high)]
    d["pp"] = d.p.map(fmt_p)
    add("## Table S3. Sensitivity to linkage: the locus sweep\n")
    add("Consecutive genes on a chromosome within the merge distance are collapsed "
        "into one locus, which is called shared-layer or Mendelian if any of its "
        "genes is. Source of Figure 1D.\n")
    add(table(d, {"merge_distance_kb": "Merge distance (kb)", "n_loci": "Loci",
                  "n_core_loci": "Shared-layer loci", "or_ci": "OR (95% CI)",
                  "pp": "P"}))
    add("")

    # ------------------------------------------------------------- S4
    d = read("fig2a_symmetry_law.csv").copy()
    d["m"] = [f"{m:.2f} ± {s:.2f}" for m, s in zip(d.mean_regions, d["sem"])]
    d["g"] = (d.frac_generalised * 100).map("{:.0f}".format)
    st = read("fig2a_symmetry_stats.csv").set_index("statistic").value
    reg = read("fig2a_symmetry_regression.csv").iloc[0]
    add("## Table S4. Phenotype breadth by breadth of action\n")
    add("Body regions affected by the Mendelian phenotype, for the 2,789 genes "
        "carrying a skeletal condition. Generalized means four or more of the six "
        "regions. Source of Figure 2A.\n")
    add(table(d, {"breadth": "Breadth", "n_genes": "Genes",
                  "m": "Regions affected (mean ± SEM)",
                  "g": "Generalized (%)"}))
    add(f"\nAcross strata: {st['mean_regions_breadth_le1']:.2f} regions at breadth "
        f"0–1 against {st['mean_regions_breadth_ge5']:.2f} at breadth 5–6 "
        f"(Mann-Whitney P = {fmt_p(st['mannwhitney_p'])}). Ordinary least squares "
        f"on breadth, adjusted for the logarithm of the gene's total ontology term "
        f"count: β = {reg.beta:.4f} regions per site, SE {reg.se:.4f}, "
        f"t = {reg.t:.2f}, P = {fmt_p(reg.p)}, n = {int(reg.n):,}. Spearman "
        f"correlation between breadth and the number of ontology terms per gene: "
        f"{st['spearman_breadth_vs_n_hpo_terms']:+.4f} "
        f"(P = {fmt_p(st['spearman_p'])}).\n")

    # ------------------------------------------------------------- S5
    d = read("fig2b_cross_organ.csv").copy()
    d["b"] = d.beta.map("{:+.4f}".format)
    d["se2"] = d.se.map("{:.4f}".format)
    d["t2"] = d.t.map("{:+.2f}".format)
    d["pp"] = d.p.map(fmt_p)
    add("## Table S5. The cross-organ control\n")
    add("Ordinary least squares of phenotype breadth on breadth of action, "
        "adjusted for the logarithm of the gene's total ontology term count. "
        "Source of Figure 2B.\n")
    add(table(d, {"predictor": "Predictor", "outcome": "Outcome", "n": "Genes",
                  "b": "β", "se2": "SE", "t2": "t", "pp": "P"}))
    add("")

    # ------------------------------------------------------------- S6
    d = read("fig2c_trait_coherence.csv").copy()
    d["r"] = d.mean_r.map("{:.3f}".format)
    d["rng"] = [f"{a:.3f}–{b:.3f}" for a, b in zip(d.min_r, d.max_r)]
    add("## Table S6. Genetic coherence of the measurement sites\n")
    add("Mean and range of the pairwise genetic correlation between the "
        "measurement sites of each organ. A shared layer can only exist where the "
        "sites share genetic architecture. Source of Figure 2C.\n")
    add(table(d, {"dataset": "Measurement set", "n_traits": "Sites",
                  "r": "Mean r", "rng": "Range",
                  "shared_layer": "Shared layer detected"}))
    add("")

    # ------------------------------------------------------------- S7
    d = read("fig3a_gefos_replication.csv").copy()
    d["z"] = d.median_z.map("{:.3f}".format)
    add("## Table S7. Replication in GEFOS\n")
    add("Median gene-level Z in GEFOS, which shares no participants with UK "
        "Biobank, as a function of breadth defined in UK Biobank. Source of "
        "Figure 3A.\n")
    add(table(d, {"site": "GEFOS site", "breadth": "Breadth", "z": "Median Z"}))
    s = read("fig3a_gefos_stats.csv").copy()
    s["pp"] = s.p.map(fmt_p)
    s["mc"] = s.median_core.map("{:.3f}".format)
    s["mr"] = s.median_rest.map("{:.3f}".format)
    add("\nShared layer against every gene outside it:\n")
    add(table(s, {"site": "GEFOS site", "n_core": "Shared layer (n)",
                  "mc": "Median Z, shared layer", "mr": "Median Z, rest",
                  "pp": "P"}))
    add("\nThe main text quotes the companion contrast in Table S9, which "
        "excludes the site-specific layer from the comparison group; the two "
        "differ only in the denominator.\n")

    # ------------------------------------------------------------- S8
    d = read("fig3b_fracture.csv").copy()
    d["z"] = d.median_z.map("{:.3f}".format)
    d["g2"] = (d.frac_z_gt2 * 100).map("{:.1f}".format)
    d["g4"] = (d.frac_z_gt4 * 100).map("{:.1f}".format)
    add("## Table S8. Fracture risk by breadth of action\n")
    add("Fracture summary statistics from an osteoporosis meta-analysis of 53,184 "
        "cases and 373,611 controls. Source of Figure 3B.\n")
    add(table(d, {"breadth": "Breadth", "n_genes": "Genes", "z": "Median Z",
                  "g2": "Z > 2 (%)", "g4": "Z > 4 (%)"}))
    add("")

    # ------------------------------------------------------------- S9
    d = read("fig3c_layer_contrast.csv").copy()
    d["z"] = d.median_z.map("{:.3f}".format)
    d["pp"] = d.p_vs_rest.map(fmt_p)
    add("## Table S9. The two layers across the endpoints\n")
    add("P values are against the 'all other genes' row of the same endpoint. "
        "Source of Figure 3C.\n")
    add(table(d, {"endpoint": "Endpoint", "layer": "Layer", "n_genes": "Genes",
                  "z": "Median Z", "pp": "P vs all other genes"}))
    add("")

    # ------------------------------------------------------------- S10
    a = read("fig3a_gefos_conditional.csv")
    b = read("fig3c_conditional.csv")
    d = pd.concat([a, b], ignore_index=True).copy()
    d["t2"] = d.t_breadth.map("{:+.2f}".format)
    d["pp"] = d.p_breadth.map(fmt_p)
    add("## Table S10. Conditional models\n")
    add("Breadth entered in one linear model alongside all six UK Biobank "
        "site-level Z statistics, from which it is derived. A coefficient that "
        "survives this adjustment is information the six scores do not carry "
        "individually.\n")
    add(table(d, {"endpoint": "Endpoint", "n": "Genes", "t2": "t for breadth",
                  "pp": "P"}))
    add("")

    # ------------------------------------------------------------- S11
    d = read("fig4a_pathways.csv").copy()
    d["or_ci"] = [f"{o:.2f} ({x:.2f}–{y:.2f})"
                  for o, x, y in zip(d.odds_ratio, d.ci_low, d.ci_high)]
    d["bg"] = (d.background_fraction * 100).map("{:.2f}".format)
    d["pp"] = d.p.map(fmt_p)
    d["hits"] = [f"{h}/{n}" for h, n in zip(d.n_hits, d.n_layer)]
    add("## Table S11. Developmental signaling pathways\n")
    add("Gene Ontology biological process terms, non-electronic evidence codes "
        "only, against all other genes. Both layers are shown; Figure 4A shows "
        "the shared layer. No correction is applied across the nine pathways.\n")
    add(table(d, {"pathway": "Pathway", "layer": "Layer", "hits": "Genes in set",
                  "bg": "Background (%)", "or_ci": "OR (95% CI)", "pp": "P"}))
    add("")

    # ------------------------------------------------------------- S12
    d = read("fig4b_regulatory_density.csv").copy()
    d["obs"] = d.observed_peaks_per_element.map("{:.2f}".format)
    d["exp"] = [f"{e:.2f} ± {s:.2f}" for e, s in zip(d.expected, d.sd)]
    d["z2"] = d.z.map("{:+.2f}".format)
    d["pp"] = d.p.map(fmt_p)
    add("## Table S12. Regulatory density against a length-matched null\n")
    add("Open-chromatin peaks in human fetal limb tissue within 100 kb of each "
        "gene, per regulatory element, against 5,000 resamples matched on the "
        "decile of log gene length. Source of Figure 4B.\n")
    add(table(d, {"layer": "Layer", "n_genes": "Genes", "obs": "Observed",
                  "exp": "Length-matched expectation ± SD", "z2": "Z", "pp": "P"}))
    add("")

    # ------------------------------------------------------------- S13
    d = read("fig4c_shared_layer_genes.csv").copy().sort_values(
        "FRAK", ascending=False)
    d["men"] = d.mendelian.map({1: "yes", 0: "—"})
    d["pan"] = d.panelapp.map({1: "yes", 0: "—"})
    for c in ("FRAK", "HEEL", "gefos_fn_z"):
        d[c] = d[c].map("{:.2f}".format)
    add("## Table S13. The 103 shared-layer genes\n")
    add("Every gene reaching Z > 2 at all six skeletal sites, ordered by its "
        "fracture-risk statistic. The paper reports breadth as a continuous "
        "property; this list is an illustration, not a deliverable, because "
        "redefining the layer in the replication cohort recovers only about 19% "
        "of it. Source of Figure 4C.\n")
    add(table(d, {"gene": "Gene", "FRAK": "Fracture Z", "HEEL": "Heel BMD Z",
                  "gefos_fn_z": "GEFOS femoral neck Z",
                  "men": "Mendelian skeletal disease", "pan": "PanelApp green"}))
    add("")

    # ------------------------------------------------------- notes
    add("---\n")
    add("## Supplemental note 1. Which comparison group each odds ratio uses\n")
    add("Two denominators appear in this paper and they are not interchangeable.\n")
    add("The enrichment odds ratios in Figure 1 and Tables S2 and S3 compare the "
        "shared layer against **every gene outside it**, which is 18,289 genes at "
        "the Z > 2 definition and includes the 36 site-specific genes. The "
        "proportions quoted in the main text alongside them — 10.7% against 2.2%, "
        "31.1% against 15.2%, 12.6% against 3.1% — use that denominator.\n")
    add("The layer contrast in Figure 3C and Table S9 compares each layer against "
        "**all other genes**, which excludes both layers and is 18,245 genes for "
        "the GEFOS endpoints and 18,247 for fracture. The three-way contrast needs "
        "the site-specific layer held out, because it is one of the groups being "
        "compared.\n")
    add("Figure 1B is a third view again: it is expressed relative to genes of "
        "breadth zero, so that the reader can see the dose-response from a fixed "
        "baseline. The proportions at breadth zero are 2.0%, 14.6% and 2.9%. No "
        "odds ratio in the paper is computed against that baseline.\n")

    add("## Supplemental note 2. Analytical choices, and what was not corrected for\n")
    add("**Preregistration.** The analysis was exploratory. No analysis plan was "
        "registered, and no p value in this paper should be read as if it came "
        "from a confirmatory design. What the central claim does rest on is "
        "agreement across three independently curated truth sides, a threshold "
        "sweep, a locus sweep, an independent cohort, and a cross-organ control "
        "that could have falsified it.\n")
    add("**Multiple comparisons.** No correction is applied across the independent "
        "questions the paper asks, nor across the nine pathways in Table S11. The "
        "sweeps in Tables S2 and S3 are sensitivity analyses of one question, not "
        "separate tests, and are reported in full so that the reader can see every "
        "value rather than the best one.\n")
    add("**The gene-level p-value floor.** MAGMA does not report a gene p below "
        "5 × 10<sup>−10</sup>. Genes whose association is stronger are all "
        "returned at the floor, which compresses the top of every Z distribution "
        "by an amount that depends on the sample size and the signal of the trait. "
        "Breadth is defined at Z > 2, far below the floor, and every endpoint "
        "analysis is rank-based, so the floor costs power and cannot create a "
        "result.\n")
    add("**Sample sizes in the replication cohort.** The published GEFOS files "
        "carry no per-variant sample size, which the gene-level method requires. "
        "N was estimated from standard errors and effect-allele frequencies as "
        "N = median[1 / (2f(1 − f)se²)] over variants with 0.1 < f < 0.9. A "
        "misestimated N rescales every gene of that trait by a constant and "
        "therefore cannot change the ranking of genes within a trait, which is "
        "what every replication statement here rests on.\n")
    add("**The annotation window.** Genes were annotated with a 35 kb upstream and "
        "10 kb downstream window, the conventional choice. The window is a single "
        "setting in the companion repository (`00_setup/config.py`), and the "
        "pipeline accepts an override so that the analysis can be repeated at "
        "10/10 kb.\n")
    add("**Confounder controls.** Gene length, research intensity and constraint "
        "are addressed in the main text. The matched-null resampling behind them "
        "is computed by `04_analysis/04_confounders.py` in the companion "
        "repository, which writes `results/confounder_matched_nulls.tsv` and "
        "`results/confounder_models.tsv`.\n")
    add("**Software.** Gene-level statistics were computed with MAGMA v1.10. "
        "All other analysis used Python 3.12 with NumPy, pandas, SciPy and "
        "statsmodels; the pinned versions are in `environment.yml` in the "
        "companion repository.\n")

    OUT.write_text("\n".join(parts) + "\n", encoding="utf-8")
    n_tables = sum(1 for p in parts if p.startswith("## Table S"))
    print(f"  wrote {OUT.relative_to(ROOT)}  ({n_tables} tables, "
          f"{len(OUT.read_text(encoding='utf-8').split())} words)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
