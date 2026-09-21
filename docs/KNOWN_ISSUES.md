# Known issues

The limitations of this analysis, stated by us. Each entry says what the
problem is, how large it is, where in the code it shows up, and what we did
about it. None of them is hidden in a supplement; the scripts print the
relevant ones as they run, and the ones that bear on a published claim are
repeated in the paper's own supplemental information.

Every figure in this file was re-checked against the tables in `results/` on
21 September 2026. Where a value here and a value in `results/` disagree, the
table is right and this file is the bug.

---

## 1 · MAGMA does not report a gene p below 5e-10

MAGMA's gene-level p-value has a numerical floor at 5e-10. Genes whose true
association is stronger are all reported at the floor, which compresses the top
of every Z distribution. The effect is not uniform across traits: it scales
with sample size and signal. The share of genes sitting on the floor is
computed per trait rather than assumed, and is heaviest for the largest
studies.

*Consequence.* Any statistic that reads the top of the distribution — a mean, a
variance, a least-squares fit on raw Z — is biased where the truncation is
heavy. We found this after computing medians on the height trait that made no
sense.

*What we do.* Every regression whose target could be affected is run on
**percentile ranks**, not raw Z: see `04_analysis/05_replication_gefos.py` and
`04_analysis/06_clinical_endpoints.py`. The share of genes on the floor is
computed per trait and written to `results/gene_coverage.tsv` by
`01_gene_scores/03_assemble_gene_matrix.py`, so the exposure of each trait is
visible rather than assumed.

*Residual risk.* Breadth itself is defined at Z > 2, far below the floor, and
is therefore unaffected. The endpoint analyses are rank-based. What remains is
that a genuinely enormous effect cannot be distinguished from a merely large
one, which costs power and cannot inflate a result.

---

## 2 · The GEFOS sample sizes are estimated, not read

The published GEFOS/UK10K files carry no per-variant or per-trait sample size,
and MAGMA needs one. N is estimated from the standard errors and effect-allele
frequencies of common variants:

    N = median[ 1 / (2 f (1-f) se^2) ],  over 0.1 < f < 0.9

giving 34,806 (femoral neck), 25,759 (lumbar spine) and 8,112 (forearm). The
estimator is standard and reproduces the published cohort sizes to within a few
percent where they are stated.

*Consequence.* MAGMA's gene Z scales with sqrt(N), so a misestimated N shifts
every gene of that trait by a constant factor. That changes absolute Z values;
it cannot change the *ranking* of genes within a trait, and every replication
statement we make is about ranks and about contrasts between gene groups within
the same trait.

*What we do.* The estimate is computed in
`01_gene_scores/01_harmonise_sumstats.py`, printed, and written to
`sample_sizes.tsv` with its provenance recorded as "estimated from se and EAF"
rather than "published".

---

## 3 · Redefining the structure in GEFOS recovers only about 19 % of the genes

This is the most important limitation in the paper and it shaped how the claim
is worded.

The *property* replicates strongly: the shared layer defined in UK Biobank
shows median Z 2.43 and 2.12 in GEFOS against 0.11 and 0.08 for all other
genes (p 4.5e-39 and 3.1e-32, `fig3c_layer_contrast.csv`), and breadth holds
at t +11.52 and +12.42 with all six UK Biobank site Z values in the model
(`fig3a_gefos_conditional.csv`).

But if the two layers are defined **afresh inside GEFOS**, the resulting
shared layer overlaps the UK Biobank one by only about **19 %**, and the
Mendelian enrichment of the re-derived layer is weaker. That enrichment is
printed by `05_replication_gefos.py` and written to
`results/replication_redefined.tsv`; it is not among the panel tables and is
therefore not quoted here.

*Consequence.* A specific list of "core genes" is not a stable object. Which
genes cross a threshold at six sites depends on the power of each site in each
cohort, and those differ.

*What we do.* The paper reports breadth as a **continuous property** and never
as a gene list. The core and the site-specific layer exist in the code because
a contrast needs sides, and the core is shown as an illustration; no claim
rests on its membership. The analysis lives in
`04_analysis/05_replication_gefos.py`, which prints the overlap and writes
`results/replication_overlap.tsv`.

---

## 4 · At locus level the direct layer contrast on fracture loses significance

Collapsing genes into loci is the right way to answer "is this just linkage
disequilibrium". The enrichment survives it: 103 shared-layer genes become 98
loci at 0 bp and 42 at 1 Mb, and the odds ratio against background *rises*
with merge distance rather than collapsing — 2.45 at the gene level to 29.5 at
1 Mb, p < 1e-4 throughout (`fig1d_locus_sweep.csv`). The rise is expected:
merging shrinks the background faster than the layer, and the confidence
intervals widen with it, from 1.6-3.7 to 1.8-485.

What does not survive is the **direct contrast between the two layers** on
fracture: at locus level it reaches only **p ≈ 0.20**.

*Consequence.* The dissociation between the shared and the site-specific layer
on fracture is a gene-level statement. At locus level it is directionally the
same and not significant.

*Why.* The site-specific layer is small to begin with (36 genes), and merging
leaves too few independent loci to contrast against. This is a power limit, not
a contradiction — but it is a power limit we cannot argue away, and the paper
states the p value rather than only the enrichment that survives.

*Where.* `02_layers/02_locus_clumping.py` computes it, prints the warning, and
writes `results/locus_endpoint_contrast.tsv`. The full sweep is in
`results/enrichment_locus_sweep.tsv`.

---

## 5 · The result depends on the MAGMA annotation window

Genes are annotated with a 35 kb upstream, 10 kb downstream window — the
conventional choice, which captures most proximal regulatory variation. It is a
choice, and a different window produces a different gene set, because a wider
window pulls more variants into long genes and into gene-dense regions.

*What we do.* The window is a single setting in `00_setup/config.py`
(`MAGMA_WINDOW`), and `01_gene_scores/02_run_magma.sh` accepts an override so
that the whole pipeline can be re-run at 10/10 kb:

    MAGMA_WINDOW=10,10 bash 01_gene_scores/02_run_magma.sh --suffix _w10

The sensitivity analysis is reported in the supplement. The direction and
significance of every headline result are unchanged; the absolute size of the
core is not, which is the same point as issue 3.

*Residual risk.* We have tested two windows, not a continuum. A pathological
window could in principle behave differently, and we have no argument that
excludes it.

---

## 6 · Binary against quantitative: the endpoint pattern is suggestive, not shown

*Not reported in the paper.* The endpoint comparison below is computed by the
pipeline but was cut from the manuscript; it is kept here because the pipeline
still produces it and a reader of the code will meet it.

Across five endpoints, breadth leads at the two **disease** endpoints (fracture
t +7.22 against peak +6.49; otosclerosis +4.66 against +2.65) and peak strength
leads at the three **measurements** (heel, height, sitting-height ratio). A
clean two-against-three split, in the expected direction.

*The alternative explanation we cannot exclude.* The two disease endpoints are
binary and the three measurements quantitative. Binary traits have lower
effective sample size for the same N, which systematically shifts the
peak-to-breadth ratio in exactly the observed direction — with no biology
involved. And a perfect 2-against-3 split occurs by chance with **p = 0.10**.

*What we do.* The paper reports the pattern as suggestive and names this
alternative in the same paragraph, and
`04_analysis/06_clinical_endpoints.py` prints it next to the table so that
nobody reads the result without it. Figure 3C carries the caveat in its legend.

---

## 7 · Otosclerosis does not separate the two layers

*Not reported in the paper.* Otosclerosis rests on a non-public input and no
longer appears in the manuscript.

Breadth predicts otosclerosis (core 15.5 % above Z 2 against 5.5 % background,
p 2.6e-05; breadth t +4.42 with all six site Z in the model). But the **layer
contrast fails**: core against site-specific gives p 0.93, and at the extreme
tail the site-specific layer carries *more* signal (8.3 % above Z 4 against
2.9 %).

The same is true of sitting-height ratio, where the layers do not separate
(p 0.90).

*What we do.* Both are reported. An endpoint that supports breadth but fails
the layer contrast is evidence about breadth and evidence against an
over-strong reading of the two-layer picture, and the paper uses it as both.

---

## 8 · Same layer, not the same genes

Within the shared layer, the fracture signal and Mendelian disease are **not**
associated (Fisher OR 1.54, p 0.40). Twenty of the 103 genes do both,
thirty-five carry only a fracture signal, thirteen only Mendelian disease, and
thirty-five carry neither (`fig4c_shared_layer_genes.csv`).

*Consequence.* The framing "one gene property, two clinical ends" is about a
**layer**, not about individual genes. The sentence "the same gene that causes
the child's dysplasia causes the grandmother's fracture" would be false and we
do not write it. What is true is that both burdens sit in the shared layer and
neither sits in the site-specific one.

---

## 9 · This analysis was exploratory, without pre-registration

Thirty-seven tests were computed in the course of this work. Twelve are
reported. Without a pre-registered plan, that ratio is part of how every p
value here should be read, and the supplement of the paper states it as a
paragraph of text.

*What that does and does not mean.* It means no single p value in this paper
should be treated as if it came from a confirmatory design. It does not mean
the results are a selection of noise: the central claims survive three
independent curations, a threshold sweep, a locus sweep, four matched-null
confounder controls, an independent cohort, and a cross-organ control that
could have killed them and did not.

*What it does not excuse.* Breadth was defined once, on the six DXA sites, and
never adjusted to improve an outcome. The layer thresholds in
`00_setup/config.py` are the ones used from the first analysis onward, and the
sweep in `02_layers/03_threshold_sweep.py` exists so that a reader can check
that nothing hinges on them.

---

## 10 · Scope

The result is a statement about the human skeleton. Three attempts to generalise
it to other organ systems are blocked for structural reasons — the subcortical
volumes share too little signal between regions to have a shared layer at all
(mean inter-region correlation 0.12 against 0.44 for the skeleton,
`fig2c_trait_coherence.csv`), and
matching GWAS anatomy to ontology anatomy fails elsewhere. The cross-organ
control in `04_analysis/03_cross_organ_control.py` reports that correlation
explicitly, because "the law does not hold in the brain" and "the brain has no
shared layer for it to hold in" are different statements and only the second is
supported.

---

## 11 · The analysis stage does not regenerate the figure inputs

The figures and the manuscript self-test read the panel tables
`results/fig*.csv`. Nothing in this repository writes them. `make analysis`
writes its own tables — `enrichment_*.tsv`, `symmetry_*.tsv`,
`cross_organ_*.tsv`, `replication_*.tsv`, `endpoints_*.tsv`,
`pathway_enrichment.tsv`, `regulatory_density.tsv` — and the `fig*.csv` files
that ship here are the frozen outputs of the run the paper reports.

*Consequence.* `make figures` and `make verify` work from a clean checkout, and
every number in the paper can be traced to a file. But re-running the analysis
does not reproduce those files, so a reader cannot yet confirm that the frozen
panel tables follow from the pipeline. This is the largest reproducibility gap
in the repository, and nothing about the figures makes it visible.

*Why it is not a rename.* The panel tables carry quantities the analysis tables
do not:

| panel table | needs, beyond the analysis table |
|---|---|
| `fig1b_dose_response.csv` | per-cell positive counts and 95% Jeffreys intervals |
| `fig1b_trend_tests.csv` | a logistic regression of the curation indicator on breadth; the analysis stage computes a Spearman correlation instead |
| `fig1c_threshold_sweep.csv` | Wald intervals on the Haldane odds ratios |
| `fig1d_locus_sweep.csv` | Wald intervals on the Haldane odds ratios |
| `fig2a_symmetry_stats.csv` | the Mann-Whitney comparison of the extreme strata, as one row |
| `fig3c_layer_contrast.csv` | each layer against *all other genes*, excluding the other layer, rather than against every non-layer gene |

*What has to be written.* One script, `04_analysis/08_export_panels.py`, that
reads the analysis tables, computes the six quantities above from the counts
they already contain, and writes every `fig*.csv`. It must be the only writer
of those files, and `make figures` must depend on it. Adding it would close the
loop from raw summary statistics to the manuscript.

---

## 12 · The limb ATAC accession was wrong; Figure 4B must be regenerated

`00_setup/02_download_references.sh` fetched **GSE170199** and described it as
fetal chromatin accessibility. GSE170199 is a two-sample transcription-factor
ChIP-seq series in HepG2 from ENCODE. Anyone who ran `make data` obtained peak
files of the wrong assay, and the regulatory-density measure behind Figure 4B
therefore has no recorded provenance here.

*What the dataset actually is.* **GSE252289**, the ATAC-seq series of Richard
et al., *Cell* 2025;188:15-32.e24. It is the only published resource matching
the description in `07_pathways_and_axis.py` — human embryonic limb cartilage,
hg19, and eight skeletal elements, namely the proximal and distal end of the
femur, tibia, humerus and radius. `00_setup/04_prepare_atac_peaks.py` now
downloads the published peak workbook and writes the eight elements as the
BED files the analysis reads; the existing loader needed no change.

*What a partial recomputation shows.*
`04_analysis/07b_recheck_regulatory_density.py` redoes the shared-layer half
of the panel directly from GSE252289, without the pipeline, and puts the two
side by side:

| | observed | length-matched null | Z |
|---|---:|---:|---:|
| recomputed from GSE252289 | 7.50 | 5.50 ± 0.46 | +4.39 |
| shipped in `fig4b_regulatory_density.csv` | 13.16 | 11.11 ± 0.55 | +3.69 |

The **absolute densities are not reproduced** — they are roughly 1.8-fold
higher in the shipped table, which is what one expects from a different peak
set. The **effect is**: an excess of about two peaks per element over a
length-matched null, in the same direction, at a similar Z. The relative
excess is in fact larger with the correct data, 1.36-fold against 1.18-fold.
The conclusion the panel supports therefore survives the correction; the
numbers printed on it do not.

*Residual risk, and it is the open one.* Identifying the right dataset does
not retroactively make the shipped numbers come from it, and the recheck
script cannot close the gap on its own: it has to use every gene in the MAGMA
location table as the background instead of the 18,392 genes that carried a Z
at all six sites, and the 36 site-specific genes — the second point of the
panel — are not recoverable from what ships here at all. **Figure 4B has to be
regenerated from the full pipeline against GSE252289 before the paper is
submitted.** Until then the observed value, the null and the Z printed in that
panel are unverified.

*Also worth fixing while there.* The `p` column of
`fig4b_regulatory_density.csv` is a two-sided normal p derived from Z, not the
empirical p over the resamples that `07_pathways_and_axis.py` writes, and the
Methods say 5,000 resamples where the script's default is 2,000
(`N_PERMUTATIONS` in `00_setup/config.py`). Both are further signs that the
panel table did not come from the pipeline as it now stands; see issue 11.
