# Supplemental information

**How widely a gene acts across the skeleton predicts how widely its mutations cause disease**

Alexander Schulz and Christian T. Thiel

Institute of Human Genetics, Universitätsklinikum Erlangen, Friedrich-Alexander-Universität Erlangen-Nürnberg, 91054 Erlangen, Germany  
Correspondence: Christian.Thiel@uk-erlangen.de

Every table below is generated from the result tables in `results/` by `06_manuscript/build_supplement.py`, so none of it can disagree with a figure or with the main text.

---

## Contents

Table S1. Breadth against Mendelian skeletal disease, three curations, breadth 0 to 6.  
Table S2. Sensitivity of the enrichment to the Z threshold that defines action at a site.  
Table S3. Sensitivity of the enrichment to collapsing genes into loci.  
Table S4. Phenotype breadth by breadth of action.  
Table S5. The cross-organ control.  
Table S6. Genetic coherence of the measurement sites within each organ.  
Table S7. Replication in GEFOS.  
Table S8. Fracture risk by breadth of action.  
Table S9. The two layers across replication and fracture endpoints.  
Table S10. Conditional models: breadth beside all six site scores.  
Table S11. Developmental signaling pathways in both layers.  
Table S12. Regulatory density against a length-matched null.  
Table S13. The 103 shared-layer genes.  
Supplemental note 1. Which comparison group each odds ratio uses.  
Supplemental note 2. Analytical choices, and what was not corrected for.  

---

## Table S1. Breadth against Mendelian skeletal disease

Proportion of genes carrying an annotation in each curation, by breadth of action. Intervals are 95% Jeffreys intervals. This is the source of Figure 1B.

| Curation | Breadth | Genes | Annotated | Annotated (%) | 95% CI (%) |
|---|---|---|---|---|---|
| PanelApp | 0 | 14242 | 290 | 2.0 | 1.8–2.3 |
| PanelApp | 1 | 2284 | 63 | 2.8 | 2.1–3.5 |
| PanelApp | 2 | 833 | 16 | 1.9 | 1.1–3.0 |
| PanelApp | 3 | 460 | 15 | 3.3 | 1.9–5.2 |
| PanelApp | 4 | 264 | 8 | 3.0 | 1.4–5.6 |
| PanelApp | 5 | 206 | 13 | 6.3 | 3.6–10.3 |
| PanelApp | 6 | 103 | 11 | 10.7 | 5.8–17.7 |
| HPO | 0 | 14242 | 2075 | 14.6 | 14.0–15.2 |
| HPO | 1 | 2284 | 390 | 17.1 | 15.6–18.7 |
| HPO | 2 | 833 | 141 | 16.9 | 14.5–19.6 |
| HPO | 3 | 460 | 84 | 18.3 | 14.9–22.0 |
| HPO | 4 | 264 | 48 | 18.2 | 13.9–23.2 |
| HPO | 5 | 206 | 43 | 20.9 | 15.8–26.8 |
| HPO | 6 | 103 | 32 | 31.1 | 22.7–40.4 |
| ClinVar | 0 | 14242 | 409 | 2.9 | 2.6–3.2 |
| ClinVar | 1 | 2284 | 87 | 3.8 | 3.1–4.7 |
| ClinVar | 2 | 833 | 31 | 3.7 | 2.6–5.2 |
| ClinVar | 3 | 460 | 21 | 4.6 | 2.9–6.8 |
| ClinVar | 4 | 264 | 11 | 4.2 | 2.2–7.1 |
| ClinVar | 5 | 206 | 12 | 5.8 | 3.2–9.7 |
| ClinVar | 6 | 103 | 13 | 12.6 | 7.3–20.0 |

## Table S2. Sensitivity to the Z threshold

The shared layer redefined at each threshold, against every gene outside it. Haldane-corrected odds ratios with Wald intervals. Source of Figure 1C.

| Z threshold | Curation | Shared layer (n) | OR (95% CI) | P |
|---|---|---|---|---|
| 2.0 | PanelApp | 103 | 5.48 (2.95–10.20) | 2.1 × 10<sup>-5</sup> |
| 2.0 | HPO | 103 | 2.53 (1.67–3.84) | 4.1 × 10<sup>-5</sup> |
| 2.0 | ClinVar | 103 | 4.62 (2.59–8.24) | 2.3 × 10<sup>-5</sup> |
| 2.5 | PanelApp | 51 | 7.38 (3.38–16.10) | 1.4 × 10<sup>-4</sup> |
| 2.5 | HPO | 51 | 3.34 (1.90–5.86) | 1.0 × 10<sup>-4</sup> |
| 2.5 | ClinVar | 51 | 5.18 (2.38–11.29) | 0.001 |
| 3.0 | PanelApp | 33 | 6.64 (2.45–17.98) | 0.006 |
| 3.0 | HPO | 33 | 4.13 (2.09–8.17) | 1.7 × 10<sup>-4</sup> |
| 3.0 | ClinVar | 33 | 5.92 (2.37–14.81) | 0.004 |
| 3.5 | PanelApp | 18 | 13.51 (4.67–39.09) | 6.1 × 10<sup>-4</sup> |
| 3.5 | HPO | 18 | 8.52 (3.39–21.40) | 1.2 × 10<sup>-5</sup> |
| 3.5 | ClinVar | 18 | 12.51 (4.63–33.84) | 1.9 × 10<sup>-4</sup> |
| 4.0 | PanelApp | 14 | 8.67 (2.22–33.82) | 0.039 |
| 4.0 | HPO | 14 | 9.59 (3.35–27.44) | 4.3 × 10<sup>-5</sup> |
| 4.0 | ClinVar | 14 | 9.31 (2.81–30.91) | 0.009 |

## Table S3. Sensitivity to linkage: the locus sweep

Consecutive genes on a chromosome within the merge distance are collapsed into one locus, which is called shared-layer or Mendelian if any of its genes is. Source of Figure 1D.

| Merge distance (kb) | Loci | Shared-layer loci | OR (95% CI) | P |
|---|---|---|---|---|
| 0 | 16521 | 98 | 2.45 (1.61–3.74) | 8.1 × 10<sup>-5</sup> |
| 100 | 2840 | 55 | 4.62 (2.45–8.70) | 1.5 × 10<sup>-7</sup> |
| 250 | 1321 | 51 | 8.71 (3.29–23.04) | 1.1 × 10<sup>-8</sup> |
| 500 | 688 | 44 | 15.27 (2.98–78.40) | 3.9 × 10<sup>-7</sup> |
| 1000 | 311 | 42 | 29.46 (1.79–485.22) | 1.1 × 10<sup>-5</sup> |

## Table S4. Phenotype breadth by breadth of action

Body regions affected by the Mendelian phenotype, for the 2,789 genes carrying a skeletal condition. Generalized means four or more of the six regions. Source of Figure 2A.

| Breadth | Genes | Regions affected (mean ± SEM) | Generalized (%) |
|---|---|---|---|
| 0 | 2056 | 2.71 ± 0.04 | 31 |
| 1 | 388 | 2.79 ± 0.08 | 33 |
| 2 | 140 | 2.63 ± 0.14 | 29 |
| 3 | 84 | 3.00 ± 0.18 | 38 |
| 4 | 48 | 2.83 ± 0.24 | 33 |
| 5 | 41 | 3.54 ± 0.27 | 46 |
| 6 | 32 | 3.38 ± 0.34 | 53 |

Across strata: 2.72 regions at breadth 0–1 against 3.47 at breadth 5–6 (Mann-Whitney P = 2.2 × 10<sup>-4</sup>). Ordinary least squares on breadth, adjusted for the logarithm of the gene's total ontology term count: β = 0.0776 regions per site, SE 0.0214, t = 3.63, P = 2.9 × 10<sup>-4</sup>, n = 2,789. Spearman correlation between breadth and the number of ontology terms per gene: -0.0067 (P = 0.724).

## Table S5. The cross-organ control

Ordinary least squares of phenotype breadth on breadth of action, adjusted for the logarithm of the gene's total ontology term count. Source of Figure 2B.

| Predictor | Outcome | Genes | β | SE | t | P |
|---|---|---|---|---|---|---|
| Skeletal breadth | Skeletal phenotype breadth | 2813 | +0.0707 | 0.0213 | +3.32 | 9.3 × 10<sup>-4</sup> |
| Skeletal breadth | Brain phenotype breadth | 2172 | +0.0355 | 0.0304 | +1.17 | 0.243 |
| Brain breadth | Skeletal phenotype breadth | 2813 | +0.0199 | 0.0295 | +0.67 | 0.500 |
| Brain breadth | Brain phenotype breadth | 2172 | +0.0118 | 0.0387 | +0.30 | 0.761 |

## Table S6. Genetic coherence of the measurement sites

Mean and range of the pairwise genetic correlation between the measurement sites of each organ. A shared layer can only exist where the sites share genetic architecture. Source of Figure 2C.

| Measurement set | Sites | Mean r | Range | Shared layer detected |
|---|---|---|---|---|
| Skeleton, 6 DXA sites | 6 | 0.443 | 0.251–0.636 | yes |
| Brain, 7 subcortical volumes | 7 | 0.119 | 0.041–0.208 | no |

## Table S7. Replication in GEFOS

Median gene-level Z in GEFOS, which shares no participants with UK Biobank, as a function of breadth defined in UK Biobank. Source of Figure 3A.

| GEFOS site | Breadth | Median Z |
|---|---|---|
| Femoral neck | 0 | 0.041 |
| Femoral neck | 1 | 0.257 |
| Femoral neck | 2 | 0.448 |
| Femoral neck | 3 | 0.557 |
| Femoral neck | 4 | 0.627 |
| Femoral neck | 5 | 0.761 |
| Femoral neck | 6 | 2.428 |
| Lumbar spine | 0 | 0.017 |
| Lumbar spine | 1 | 0.221 |
| Lumbar spine | 2 | 0.386 |
| Lumbar spine | 3 | 0.427 |
| Lumbar spine | 4 | 0.567 |
| Lumbar spine | 5 | 1.084 |
| Lumbar spine | 6 | 2.116 |
| Forearm | 0 | 0.077 |
| Forearm | 1 | 0.074 |
| Forearm | 2 | 0.185 |
| Forearm | 3 | 0.136 |
| Forearm | 4 | 0.232 |
| Forearm | 5 | 0.415 |
| Forearm | 6 | 0.645 |

Shared layer against every gene outside it:

| GEFOS site | Shared layer (n) | Median Z, shared layer | Median Z, rest | P |
|---|---|---|---|---|
| Femoral neck | 103 | 2.428 | 0.107 | 4.7 × 10<sup>-39</sup> |
| Lumbar spine | 103 | 2.116 | 0.083 | 3.3 × 10<sup>-32</sup> |
| Forearm | 103 | 0.645 | 0.090 | 3.2 × 10<sup>-9</sup> |

The main text quotes the companion contrast in Table S9, which excludes the site-specific layer from the comparison group; the two differ only in the denominator.

## Table S8. Fracture risk by breadth of action

Fracture summary statistics from an osteoporosis meta-analysis of 53,184 cases and 373,611 controls. Source of Figure 3B.

| Breadth | Genes | Median Z | Z > 2 (%) | Z > 4 (%) |
|---|---|---|---|---|
| 0 | 14242 | 0.221 | 5.4 | 0.1 |
| 1 | 2284 | 0.437 | 9.2 | 0.3 |
| 2 | 833 | 0.618 | 12.0 | 2.0 |
| 3 | 460 | 0.644 | 13.5 | 2.6 |
| 4 | 264 | 0.513 | 15.2 | 1.5 |
| 5 | 206 | 0.799 | 24.3 | 8.3 |
| 6 | 103 | 2.329 | 53.4 | 23.3 |

## Table S9. The two layers across the endpoints

P values are against the 'all other genes' row of the same endpoint. Source of Figure 3C.

| Endpoint | Layer | Genes | Median Z | P vs all other genes |
|---|---|---|---|---|
| GEFOS femoral neck | Shared layer | 103 | 2.428 | 4.5 × 10<sup>-39</sup> |
| GEFOS femoral neck | Site-specific layer | 36 | 0.647 | 0.040 |
| GEFOS femoral neck | All other genes | 18245 | 0.106 | — |
| GEFOS lumbar spine | Shared layer | 103 | 2.116 | 3.1 × 10<sup>-32</sup> |
| GEFOS lumbar spine | Site-specific layer | 36 | 0.113 | 0.049 |
| GEFOS lumbar spine | All other genes | 18245 | 0.083 | — |
| Fracture | Shared layer | 103 | 2.329 | 6.3 × 10<sup>-26</sup> |
| Fracture | Site-specific layer | 36 | 0.965 | 0.010 |
| Fracture | All other genes | 18247 | 0.289 | — |

## Table S10. Conditional models

Breadth entered in one linear model alongside all six UK Biobank site-level Z statistics, from which it is derived. A coefficient that survives this adjustment is information the six scores do not carry individually.

| Endpoint | Genes | t for breadth | P |
|---|---|---|---|
| GEFOS Femoral neck | 18384 | +11.52 | 1.4 × 10<sup>-30</sup> |
| GEFOS Lumbar spine | 18384 | +12.42 | 2.8 × 10<sup>-35</sup> |
| GEFOS Forearm | 18384 | +2.85 | 0.004 |
| Fracture | 18386 | +11.24 | 3.1 × 10<sup>-29</sup> |
| Heel BMD | 18383 | +19.89 | 4.1 × 10<sup>-87</sup> |

## Table S11. Developmental signaling pathways

Gene Ontology biological process terms, non-electronic evidence codes only, against all other genes. Both layers are shown; Figure 4A shows the shared layer. No correction is applied across the nine pathways.

| Pathway | Layer | Genes in set | Background (%) | OR (95% CI) | P |
|---|---|---|---|---|---|
| FGF signaling | Shared layer | 6/103 | 0.42 | 16.75 (7.33–38.29) | 4.9 × 10<sup>-6</sup> |
| FGF signaling | Site-specific layer | 0/36 | 0.42 | 3.19 (0.19–52.43) | 1.000 |
| Ossification | Shared layer | 12/103 | 2.15 | 6.38 (3.50–11.62) | 2.1 × 10<sup>-6</sup> |
| Ossification | Site-specific layer | 3/36 | 2.15 | 4.78 (1.58–14.45) | 0.042 |
| Wnt, canonical | Shared layer | 9/103 | 1.61 | 6.29 (3.20–12.39) | 4.3 × 10<sup>-5</sup> |
| Wnt, canonical | Site-specific layer | 2/36 | 1.61 | 4.44 (1.22–16.13) | 0.114 |
| Wnt, all | Shared layer | 11/103 | 2.31 | 5.37 (2.89–10.00) | 2.6 × 10<sup>-5</sup> |
| Wnt, all | Site-specific layer | 2/36 | 2.31 | 3.08 (0.85–11.14) | 0.201 |
| Cartilage development | Shared layer | 5/103 | 1.04 | 5.39 (2.26–12.88) | 0.005 |
| Cartilage development | Site-specific layer | 2/36 | 1.04 | 6.91 (1.90–25.15) | 0.054 |
| BMP signaling | Shared layer | 1/103 | 0.78 | 1.86 (0.37–9.42) | 0.553 |
| BMP signaling | Site-specific layer | 1/36 | 0.78 | 5.40 (1.05–27.90) | 0.245 |
| Hedgehog signaling | Shared layer | 1/103 | 0.70 | 2.07 (0.41–10.47) | 0.517 |
| Hedgehog signaling | Site-specific layer | 1/36 | 0.70 | 5.99 (1.16–30.99) | 0.224 |
| TGF-beta signaling | Shared layer | 0/103 | 1.06 | 0.45 (0.03–7.22) | 1.000 |
| TGF-beta signaling | Site-specific layer | 1/36 | 1.06 | 3.95 (0.77–20.34) | 0.319 |
| Notch signaling | Shared layer | 3/103 | 0.94 | 3.72 (1.27–10.92) | 0.072 |
| Notch signaling | Site-specific layer | 2/36 | 0.94 | 7.73 (2.12–28.15) | 0.045 |

## Table S12. Regulatory density against a length-matched null

Open-chromatin peaks in human fetal limb tissue within 100 kb of each gene, per regulatory element, against 5,000 resamples matched on the decile of log gene length. Source of Figure 4B.

| Layer | Genes | Observed | Length-matched expectation ± SD | Z | P |
|---|---|---|---|---|---|
| Shared layer | 103 | 13.16 | 11.11 ± 0.55 | +3.69 | 2.2 × 10<sup>-4</sup> |
| Site-specific layer | 36 | 11.58 | 11.35 ± 0.96 | +0.24 | 0.807 |

## Table S13. The 103 shared-layer genes

Every gene reaching Z > 2 at all six skeletal sites, ordered by its fracture-risk statistic. The paper reports breadth as a continuous property; this list is an illustration, not a deliverable, because redefining the layer in the replication cohort recovers only about 19% of it. Source of Figure 4C.

| Gene | Fracture Z | Heel BMD Z | GEFOS femoral neck Z | Mendelian skeletal disease | PanelApp green |
|---|---|---|---|---|---|
| WNT16 | 7.96 | 7.95 | 3.81 | — | — |
| SLC26A1 | 7.85 | 8.03 | 4.37 | yes | — |
| FGFRL1 | 7.84 | 7.62 | 4.70 | yes | — |
| FAM3C | 7.64 | 10.21 | 3.22 | — | — |
| CPED1 | 7.21 | 10.32 | 2.27 | — | — |
| IDUA | 6.79 | 7.56 | 3.47 | yes | yes |
| CCDC170 | 6.71 | 10.12 | 7.21 | — | — |
| DGKQ | 6.64 | 7.01 | 3.22 | — | — |
| PPP6R3 | 6.03 | 7.20 | 4.95 | — | — |
| LRP5 | 5.20 | 9.45 | 4.25 | yes | yes |
| SOX6 | 5.04 | 10.63 | 5.54 | yes | — |
| SHFM1 | 4.92 | 6.11 | 5.60 | — | — |
| F2 | 4.89 | 6.11 | 4.10 | — | — |
| PACSIN3 | 4.68 | 7.20 | 3.84 | — | — |
| C11orf49 | 4.65 | 7.95 | 3.51 | — | — |
| CKAP5 | 4.59 | 7.75 | 2.68 | — | — |
| ARFGAP2 | 4.57 | 6.11 | 3.89 | — | — |
| ZNF408 | 4.57 | 7.21 | 4.45 | yes | — |
| ARHGAP1 | 4.56 | 6.11 | 4.09 | — | — |
| GAK | 4.51 | 7.28 | 1.83 | — | — |
| PIGL | 4.35 | -0.28 | -0.29 | yes | — |
| ESR1 | 4.31 | 11.27 | 4.22 | yes | — |
| C11orf58 | 4.29 | 6.11 | 3.95 | — | — |
| C7orf76 | 4.24 | 7.96 | 7.72 | — | — |
| TMEM175 | 3.95 | 7.13 | 1.38 | — | — |
| DGKZ | 3.94 | 7.37 | 3.03 | — | — |
| ATG13 | 3.83 | 6.11 | 3.70 | — | — |
| DDB2 | 3.78 | 6.11 | 3.31 | yes | — |
| MDK | 3.75 | 6.11 | 3.34 | — | — |
| WLS | 3.69 | 10.10 | 6.86 | yes | — |
| AMBRA1 | 3.63 | 8.46 | 3.76 | — | — |
| HARBI1 | 3.62 | 8.10 | 3.60 | — | — |
| PHLDB1 | 3.59 | 6.56 | 0.36 | yes | — |
| CHRM4 | 3.56 | 6.92 | 3.45 | — | — |
| IBSP | 3.45 | 9.07 | 3.52 | — | — |
| PFDN5 | 3.40 | 5.91 | 4.30 | — | — |
| CREB3L1 | 3.30 | 7.94 | 1.75 | yes | yes |
| C12orf10 | 3.22 | 5.73 | 4.08 | — | — |
| RUNX2 | 3.17 | 10.63 | 1.09 | yes | yes |
| CYP19A1 | 3.14 | 9.52 | 4.31 | yes | — |
| ESPL1 | 3.11 | 6.08 | 4.27 | — | — |
| WNT4 | 3.02 | 7.24 | 4.19 | yes | — |
| FUBP3 | 3.01 | 8.59 | 3.96 | — | — |
| JAG1 | 2.81 | 7.62 | 3.21 | yes | — |
| MEPE | 2.70 | 8.66 | 5.09 | — | — |
| CDC42 | 2.64 | 8.47 | 2.87 | yes | — |
| ASPSCR1 | 2.53 | 8.03 | 2.03 | — | — |
| NARS2 | 2.38 | 3.34 | 0.02 | yes | — |
| PYCR1 | 2.37 | 6.11 | 1.67 | yes | yes |
| NR1H3 | 2.36 | 7.73 | 2.35 | — | — |
| TNPO1 | 2.34 | 4.50 | 2.02 | — | — |
| MYADML2 | 2.33 | 6.11 | 1.73 | — | — |
| MAFG | 2.16 | 5.88 | 1.50 | — | — |
| NOTUM | 2.15 | 6.11 | 1.78 | — | — |
| CEP120 | 2.12 | 7.92 | 2.33 | yes | yes |
| PKDCC | 2.00 | 7.72 | 3.97 | yes | yes |
| TMEM171 | 1.91 | 4.24 | 1.38 | — | — |
| TNFRSF11A | 1.90 | 8.11 | 1.77 | yes | yes |
| SEMA4G | 1.85 | 3.35 | -0.38 | — | — |
| GAB2 | 1.83 | 3.18 | -0.01 | — | — |
| FCHO2 | 1.74 | 3.96 | 1.68 | — | — |
| KLHL17 | 1.74 | 6.62 | 0.62 | — | — |
| C17orf82 | 1.71 | 6.11 | 2.43 | — | — |
| SLC25A13 | 1.58 | 9.07 | 3.47 | — | — |
| C10orf2 | 1.57 | 3.41 | -0.20 | — | — |
| RAP1GAP | 1.54 | 7.94 | 0.62 | — | — |
| MICA | 1.42 | 3.88 | 1.06 | — | — |
| C2orf91 | 1.40 | 9.82 | 1.60 | — | — |
| XKR6 | 1.32 | 9.60 | 0.81 | — | — |
| TBX2 | 1.31 | 5.82 | 2.45 | yes | — |
| LOC283278 | 1.27 | 7.75 | 2.95 | — | — |
| LIN7C | 1.17 | 7.16 | 0.55 | — | — |
| CCND1 | 1.05 | 4.75 | 2.86 | yes | — |
| NFIB | 1.02 | 0.19 | 0.32 | yes | — |
| PPP1CB | 0.97 | 7.69 | 2.35 | yes | — |
| SNX24 | 0.90 | 2.75 | 2.75 | — | — |
| COLEC10 | 0.81 | 6.51 | 6.40 | yes | — |
| EIF2B2 | 0.75 | 6.99 | 2.10 | — | — |
| SPP1 | 0.68 | 7.79 | 2.06 | — | — |
| LZTS3 | 0.65 | 1.50 | -0.27 | — | — |
| CSNK1G3 | 0.55 | 7.33 | 1.50 | — | — |
| PPIC | 0.49 | 2.10 | 3.34 | — | — |
| NR1H4 | 0.43 | 0.70 | 0.13 | — | — |
| SEMA3D | 0.40 | 3.35 | 0.31 | — | — |
| TRIM42 | 0.39 | 0.41 | 1.57 | — | — |
| HIST1H2BA | 0.36 | 6.28 | 0.55 | — | — |
| TNFRSF11B | 0.27 | 7.24 | 6.62 | yes | yes |
| TANGO6 | 0.18 | 2.38 | -0.04 | — | — |
| CASR | 0.03 | 1.26 | 2.45 | yes | yes |
| MARC1 | -0.11 | 2.27 | 0.73 | — | — |
| ZNF592 | -0.20 | 3.80 | 1.61 | yes | — |
| PGF | -0.28 | 7.25 | 1.48 | — | — |
| ALPK3 | -0.33 | 3.51 | 1.58 | yes | — |
| C20orf194 | -0.39 | 2.62 | 0.57 | — | — |
| SLC4A11 | -0.42 | 2.37 | 0.22 | — | — |
| DTNBP1 | -0.69 | 2.63 | 0.08 | — | — |
| COL4A2 | -0.71 | 2.54 | 1.44 | — | — |
| POM121L2 | -0.73 | 4.14 | -0.87 | — | — |
| GALNT3 | -0.74 | 2.34 | 5.83 | yes | yes |
| SLC30A10 | -0.83 | 6.11 | 2.09 | yes | — |
| CSRNP3 | -0.90 | 3.00 | 4.94 | — | — |
| MCC | -1.22 | 2.22 | 0.53 | — | — |
| NAGK | -1.57 | 1.51 | -0.62 | — | — |

---

## Supplemental note 1. Which comparison group each odds ratio uses

Two denominators appear in this paper and they are not interchangeable.

The enrichment odds ratios in Figure 1 and Tables S2 and S3 compare the shared layer against **every gene outside it**, which is 18,289 genes at the Z > 2 definition and includes the 36 site-specific genes. The proportions quoted in the main text alongside them — 10.7% against 2.2%, 31.1% against 15.2%, 12.6% against 3.1% — use that denominator.

The layer contrast in Figure 3C and Table S9 compares each layer against **all other genes**, which excludes both layers and is 18,245 genes for the GEFOS endpoints and 18,247 for fracture. The three-way contrast needs the site-specific layer held out, because it is one of the groups being compared.

Figure 1B is a third view again: it is expressed relative to genes of breadth zero, so that the reader can see the dose-response from a fixed baseline. The proportions at breadth zero are 2.0%, 14.6% and 2.9%. No odds ratio in the paper is computed against that baseline.

## Supplemental note 2. Analytical choices, and what was not corrected for

**Preregistration.** The analysis was exploratory. No analysis plan was registered, and no p value in this paper should be read as if it came from a confirmatory design. What the central claim does rest on is agreement across three independently curated truth sides, a threshold sweep, a locus sweep, an independent cohort, and a cross-organ control that could have falsified it.

**Multiple comparisons.** No correction is applied across the independent questions the paper asks, nor across the nine pathways in Table S11. The sweeps in Tables S2 and S3 are sensitivity analyses of one question, not separate tests, and are reported in full so that the reader can see every value rather than the best one.

**The gene-level p-value floor.** MAGMA does not report a gene p below 5 × 10<sup>−10</sup>. Genes whose association is stronger are all returned at the floor, which compresses the top of every Z distribution by an amount that depends on the sample size and the signal of the trait. Breadth is defined at Z > 2, far below the floor, and every endpoint analysis is rank-based, so the floor costs power and cannot create a result.

**Sample sizes in the replication cohort.** The published GEFOS files carry no per-variant sample size, which the gene-level method requires. N was estimated from standard errors and effect-allele frequencies as N = median[1 / (2f(1 − f)se²)] over variants with 0.1 < f < 0.9. A misestimated N rescales every gene of that trait by a constant and therefore cannot change the ranking of genes within a trait, which is what every replication statement here rests on.

**The annotation window.** Genes were annotated with a 35 kb upstream and 10 kb downstream window, the conventional choice. The window is a single setting in the companion repository (`00_setup/config.py`), and the pipeline accepts an override so that the analysis can be repeated at 10/10 kb.

**Confounder controls.** Gene length, research intensity and constraint are addressed in the main text. The matched-null resampling behind them is computed by `04_analysis/04_confounders.py` in the companion repository, which writes `results/confounder_matched_nulls.tsv` and `results/confounder_models.tsv`.

**Software.** Gene-level statistics were computed with MAGMA v1.10. All other analysis used Python 3.12 with NumPy, pandas, SciPy and statsmodels; the pinned versions are in `environment.yml` in the companion repository.

