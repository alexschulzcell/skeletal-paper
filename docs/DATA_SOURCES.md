# Data sources

Every input the pipeline reads, with its sample size, accession, version,
access date and licence. Nothing in this table is distributed with the
repository; `00_setup/01_download_sumstats.sh` and
`00_setup/02_download_references.sh` fetch all of it except the two rows marked
**not public**, and `00_setup/03_check_inputs.py` verifies what arrived.

Access date for every row: **19 September 2026**, unless stated otherwise.

---

## 1 · The six sites that define breadth

Breadth is the number of these six at which a gene reaches Z > 2. They are one
DXA release from one cohort, which is the point: the six measurements are
anatomically named, measured on the same participants, on the same instrument,
in the same analysis. Any cross-cohort assembly would confound anatomy with
cohort.

| Trait | N | Accession | Source | Build | Licence |
|---|---:|---|---|---|---|
| DXA BMD, head | 31,986 | GCST90568450 | UK Biobank, via GWAS Catalog | GRCh37 | UK Biobank terms |
| DXA BMD, arms | 31,873 | GCST90568449 | UK Biobank, via GWAS Catalog | GRCh37 | UK Biobank terms |
| DXA BMD, legs | 31,873 | GCST90568451 | UK Biobank, via GWAS Catalog | GRCh37 | UK Biobank terms |
| DXA BMD, pelvis | 31,873 | GCST90568452 | UK Biobank, via GWAS Catalog | GRCh37 | UK Biobank terms |
| DXA BMD, femoral neck | 32,017 | GCST90568458 | UK Biobank, via GWAS Catalog | GRCh37 | UK Biobank terms |
| DXA BMD, lumbar spine | 30,449 | GCST90568448 | UK Biobank, via GWAS Catalog | GRCh37 | UK Biobank terms |

Source publication: Qian Y, Xia J, Wang P, et al. *Genomics Proteomics
Bioinformatics* 2025;23:qzaf097 (PMID 41206123). The sample size differs
between sites and is **not** the same number for all six; the accession and N
above were read back from the GWAS Catalog REST API on the access date. An
earlier version of this table gave 31,986 for every site and listed the
accessions in the wrong order, and the manuscript repeated the single sample
size; both are corrected.

The eleven-accession release GCST90568448–GCST90568458 contains further
sub-regions; the six above are the ones with distinct anatomical names and
non-overlapping regions of interest.

## 2 · Endpoints

These never enter the definition of breadth. Each is a target, not a predictor.

| Trait | N | Accession | Source | Modality | Licence |
|---|---:|---|---|---|---|
| Heel BMD | 426,824 | GCST006979 | Morris 2019, *Nat Genet* 51:258 | quantitative ultrasound | open, GWAS Catalog |
| Fracture | 53,184 cases, 373,611 controls | GCST006980 | Morris 2019, *Nat Genet* 51:258 | binary, disease | open, GWAS Catalog |
| Otosclerosis | 864,702 | — | **not public**, in-house harmonisation | binary, disease | see note below |
| Standing height | 1,232,747 | Yengo 2022 | GIANT consortium | quantitative | open, GIANT terms |
| Sitting-height ratio | 473,511 | GCST90728588 | UK Biobank, via GWAS Catalog | quantitative | UK Biobank terms |

**Otosclerosis.** These summary statistics are our own harmonisation of a
case/control analysis and are not deposited in a public archive. Every
otosclerosis result in the paper is marked as depending on a non-public input,
and `04_analysis/06_clinical_endpoints.py` skips that endpoint with a stated
message when the file is absent. No conclusion of the paper rests on it: it is
a disease endpoint that supports the breadth result and explicitly fails to
separate the two layers (p 0.93).

## 3 · Replication cohort

| Trait | N | Source | Overlap with discovery | Licence |
|---|---:|---|---|---|
| GEFOS femoral neck | 34,806 | Zheng 2015, GEFOS/UK10K | none | open, gefos.org |
| GEFOS lumbar spine | 25,759 | Zheng 2015, GEFOS/UK10K | none | open, gefos.org |
| GEFOS forearm | 8,112 | Zheng 2015, GEFOS/UK10K | none | open, gefos.org |

No UK Biobank participants. The published files carry no sample-size column, so
N is estimated from standard errors and allele frequencies; see
`docs/KNOWN_ISSUES.md` and `01_gene_scores/01_harmonise_sumstats.py`.

## 4 · Cross-organ control

| Trait | N | Source | Version | Licence |
|---|---:|---|---|---|
| Seven subcortical volumes (thalamus, caudate, putamen, pallidum, hippocampus, amygdala, accumbens) | 33,211 each | Oxford BIG40 | release 2, stats33k | open, BIG40 terms |

Used only in `04_analysis/03_cross_organ_control.py` and to build
`results/truth_brain.tsv`. The paper makes no claim about the brain other than
the control.

## 5 · Truth sides

Three independent curations, deliberately from unrelated pipelines: an expert
diagnostic panel, a phenotype ontology, and variant-level curation.

| Source | Content used | Version | Access | Licence |
|---|---|---|---|---|
| Genomics England PanelApp, panel 309 "Skeletal dysplasia" | green (confidence 3) genes only | as served by the API on the access date; the JSON records its own version | API, `/api/v1/panels/309/` | Apache 2.0 / open |
| Human Phenotype Ontology | `hp.obo`, `phenotype.hpoa`, `genes_to_disease.txt`; phenotype ("P") rows only | 2026-08 release | purl.obolibrary.org | HPO licence, free with attribution |
| ClinVar | `clinvar.vcf.gz`, GRCh37; ~4.5 million records; pathogenic / likely pathogenic, with `CLNDN` matched against a skeletal-diagnosis pattern; benign as negative control | weekly release, week of the access date | NCBI FTP | public domain |

The six anatomical body regions on the disease side — skull, spine, thorax,
pelvis, long bones, hands and feet — are the full HPO descendant sets of the
roots listed in `00_setup/config.py`, and are near-disjoint by construction.

## 6 · Reference files

| File | Purpose | Version | Source | Licence |
|---|---|---|---|---|
| `g1000_eur.{bed,bim,fam}` | LD reference for MAGMA | 1000 Genomes phase 3, EUR | CNCR, cncr.nl/research/magma | open |
| `NCBI37.3.gene.loc` | gene coordinates and symbols | NCBI 37.3 | CNCR | open |
| gnomAD constraint | LOEUF, for the constraint confounder | v2.1.1 | gnomAD public bucket | open, gnomAD terms |
| `gene2pubmed` | publications per gene, for the ascertainment confounder | current | NCBI FTP | public domain |
| `goa_human.gaf`, `go-basic.obo` | pathway membership; non-IEA evidence only | current | Gene Ontology Consortium | CC BY 4.0 |

## 7 · Mechanism

| Source | Content | N / scale | Licence |
|---|---|---|---|
| Human embryonic limb ATAC | open-chromatin peaks, eight limb skeletal elements, hg19 | 8 elements, E54 | GSE252289, open |
| Human growth-plate pseudobulk | zone-level expression, resting → hypertrophic | **not public** in this form; derived from a fetal skeletal atlas of 372,937 cells | see note |
| Visium, human growth plate | 16 sections, 12 donors; the tissue image of Figure 4C | 16 sections | see note |

**The limb ATAC source, and the accession that was wrong.**
`00_setup/02_download_references.sh` used to fetch **GSE170199** for this row
and describe it as "fetal chromatin accessibility". It is not: GSE170199 is a
two-sample transcription-factor ChIP-seq series in HepG2 from ENCODE — wrong
assay, wrong tissue.

The dataset the analysis describes is **GSE252289**, the ATAC-seq series of

> Richard, D., Muthuirulan, P., Young, M., Yengo, L., Vedantam, S., Marouli,
> E., Bartell, E., GIANT Consortium, Hirschhorn, J., and Capellini, T.D.
> (2025). Functional genomics of human skeletal development and the patterning
> of height heritability. *Cell* **188**, 15–32.e24.

It is the only published resource that matches the description in the analysis
exactly: human embryonic limb cartilage, peaks called on hg19, and **eight**
skeletal elements — the proximal and distal end of the femur, tibia, humerus
and radius. The workbook `GSE252289_ATAC_peaks.xlsx` carries one sheet per
element per stage; `00_setup/04_prepare_atac_peaks.py` downloads it and writes
the eight E54 sheets as `GSE252289_ATAC_<Element>_E54_hg19.bed.gz`, which is
the layout `04_analysis/07_pathways_and_axis.py` reads. Peak counts per
element run from 31,397 (ProxRadius) to 45,162 (DistFemur).

**Open point.** The numbers currently in Figure 4B were produced before this
was corrected, from peak files whose origin is not recorded here. They have to
be regenerated against GSE252289 and re-checked before the paper is submitted;
see `docs/KNOWN_ISSUES.md`, issue 12.

**The two non-public rows.** The growth-plate pseudobulk matrix and the Visium
sections are derived products of primary human tissue data. They are not
distributed here. `04_analysis/07_pathways_and_axis.py` runs its first two
sections without them and states clearly that the third is skipped; the pathway
and regulatory results, which carry the mechanism claim, do not depend on them.

---

## Reproducing the table

`results/gene_coverage.tsv`, written by
`01_gene_scores/03_assemble_gene_matrix.py`, states for every trait the sample
size actually used, the number of genes tested, and the share of genes sitting
on MAGMA's gene-p floor. If a number here and a number there disagree, the file
is right.
