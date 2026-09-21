# How widely a gene acts across the skeleton
# predicts how widely its mutations cause disease

<!-- AJHG Report. Title occupies two lines of 42 and 47 characters; the limit is three
     lines of 54 characters including spaces. Citations are superscript numerals in
     order of first appearance and are listed under References. -->

**Authors.** Alexander Schulz<sup>1</sup> (ORCID 0009-0009-2605-4350),
Christian T. Thiel<sup>1,\*</sup> (ORCID 0000-0003-3817-7277)

<sup>1</sup> Institute of Human Genetics, Universitätsklinikum Erlangen,
Friedrich-Alexander-Universität Erlangen-Nürnberg, 91054 Erlangen, Germany

<sup>\*</sup>Correspondence: Christian.Thiel@uk-erlangen.de

---

## Abstract

Mendelian disease genes crowd the loci of common-variant associations, so a gene acting on many
traits ought to cause many phenotypes when it breaks. Across the phenome it does not: the two
kinds of pleiotropy are uncorrelated. We suspected the comparison, not the biology. A
phenome-wide count treats bone density, urate and lung function as interchangeable, and weighs
them against ontology terms from every organ — which hides an anatomically local correspondence.
So we measured breadth on both sides with the same anatomy, inside one organ. Bone mineral
density at six anatomically named skeletal sites in 30,449–32,017 individuals gives each of
18,392 genes one number: how many sites it acts at. That number predicted Mendelian skeletal
disease in three independent curations — odds ratios 5.5, 2.5 and 4.6 at all six sites, and a
graded rise of 1.11 to 1.23 per site. It also predicted how much of the body the disease
involved: 3.47 of six regions against 2.72 for genes acting at one site or none, research effort
adjusted out. The control designed to kill the result did not. Rebuilt on seven brain volumes,
breadth predicted nothing in either organ: those volumes share too little genetic architecture
for a shared layer to exist — the phenome-wide null seen from inside. The layer replicated in a
cohort with no UK Biobank participants and reached the clinic, genes carrying a fracture signal
rising from 5% to 53%. One property spans two clinical ends: generalized dysplasia when severely
mutated, fracture when finely tuned.

---

## Main text

Rare and common variation in the same gene often produce related conditions, and genes carrying
Mendelian disease alleles are enriched near the association signals of common variants.<sup>1</sup>
It would follow that a gene acting on many traits in a genome-wide
association study should also cause many phenotypes when it is disrupted outright. Across the
phenome it does not. In a systematic comparison of UK Biobank association breadth against Human
Phenotype Ontology annotation breadth, the two measures were uncorrelated.<sup>2</sup>
Two firm observations therefore stand in contradiction, and the contradiction has no
accepted resolution.

We reasoned that the fault might lie in the comparison rather than in the biology. A phenome-wide
count treats bone mineral density, serum urate and forced expiratory volume as interchangeable
units of pleiotropy, and the Mendelian side counts ontology terms drawn from every organ. If the
correspondence between the two kinds of pleiotropy is anatomical rather than global, matching
counts of unrelated things will destroy it. The test that follows is to measure breadth on both
sides with the same anatomy, inside a single organ system.

The skeleton allows that test. Dual-energy X-ray absorptiometry in UK Biobank yields bone mineral
density at six sites that are named after body parts rather than after machine settings: head,
arms, legs, pelvis, femoral neck and lumbar spine, each in 30,449 to 32,017 individuals.<sup>3</sup>
The Human Phenotype Ontology names the same parts on the disease side. We computed gene-level
association statistics at each of the six sites, and gave every gene one number, its breadth: the
count of sites at which it reaches Z > 2 (Figure 1A). Of 18,392 genes, 103 act at all six sites;
we call these the shared layer. A second, smaller group of 36 genes acts strongly at exactly one
site and nowhere else, and we call these site-specific.

Breadth of action predicted Mendelian skeletal disease, and it did so in three curations that do
not share a method of ascertainment: an expert-reviewed gene panel, a phenotype ontology, and
variant-level curation (Figure 1B). Among genes acting at all six sites, 10.7% carry a
PanelApp-reviewed skeletal dysplasia gene assignment against 2.2% of the 18,289 genes that do not
(odds ratio 5.5, 95% CI 2.9–10.2, P = 2 × 10⁻⁵); 31.1% against 15.2% for a Human Phenotype
Ontology skeletal annotation (OR 2.5, 1.7–3.8, P = 4 × 10⁻⁵); and 12.6% against 3.1% for a
pathogenic ClinVar variant carrying a skeletal diagnosis (OR 4.6, 2.6–8.2, P = 2 × 10⁻⁵). Measured
against genes that act at no site at all, the same three proportions are 2.0%, 14.6% and 2.9%. The
increase is graded rather than a property of the extreme: each additional site multiplies the
odds by 1.23, 1.11 and 1.20 in the three curations (P = 9 × 10⁻⁹, 5 × 10⁻⁹ and 2 × 10⁻⁸). Benign
ClinVar variants gave no enrichment, which is the negative control the curation itself supplies.

Two artifacts could manufacture such a pattern. The first is the threshold. Defining action at a
site by Z > 2 is a choice, and we repeated the analysis at 2.5, 3.0, 3.5 and 4.0, which shrinks
the shared layer from 103 genes to 14. The odds ratio did not fall; against PanelApp it ranged
from 5.5 to 13.5 and remained significant at every cut (Figure 1C). The second is linkage. A gene
that acts at all six sites has neighbors that do the same, so 103 genes may be far fewer
independent findings. Collapsing genes into loci at merge distances from zero to one megabase
leaves 98 to 42 independent shared-layer loci, and the enrichment rises rather than falls with
distance, from 2.5 at the gene level to 4.6 at 100 kb and 15.3 at 500 kb (Figure 1D).

If a shared layer exists, it should predict not only whether a gene causes disease but how much
of the body that disease involves. We annotated 2,789 genes with a Mendelian skeletal condition
for the six body regions their phenotypes affect — skull, spine, thorax, pelvis, long bones, and
hands and feet — using the same anatomical partition as the genetic side. Genes acting at five or
six skeletal sites caused disease in 3.47 regions on average, against 2.72 for genes acting at one
site or none (P = 2 × 10⁻⁴), and the proportion whose disease is generalized, meaning four or more
regions, rose from 32% to 49% (Figure 2A).

The obvious objection is that widely studied genes accumulate more annotations of every kind. It
does not apply here. The outcome is a proportion of six anatomical regions, not a count of terms,
so research effort divides out, and the data confirm it: the correlation between breadth of action
and the total number of Human Phenotype Ontology terms per gene is −0.007. Regressing the number
of affected regions on breadth while adjusting for the logarithm of the annotation count leaves
the coefficient intact at 0.078 per site (t = 3.63, P = 3 × 10⁻⁴). Gene length gave the same
answer: the correlation between breadth and the number of variants tested per gene is −0.013, and
a length-matched null retains the enrichment. Constraint does not explain it either, since the
shared layer is not depleted of loss-of-function-tolerant genes.

The deeper objection is that breadth measures nothing skeletal — that some genes simply matter
more, everywhere, and that a gene acting at six bone sites is a gene acting on everything. This
predicts a specific result, and we tested it directly. We repeated the entire construction on
seven subcortical brain volumes in up to 33,211 individuals,<sup>4</sup> deriving a brain
breadth for every gene, and built the matching neurological phenotype breadth from the ontology.
The resulting two-by-two is unambiguous (Figure 2B). Skeletal breadth predicts skeletal phenotype
breadth (β = +0.071, t = 3.32, P = 9 × 10⁻⁴). Brain breadth predicts skeletal phenotype breadth
not at all (β = +0.020, P = 0.50). General importance would have produced a symmetric matrix; the
matrix is not symmetric.

Brain breadth also failed to predict neurological phenotype breadth (P = 0.76), and the reason
matters more than the result. A shared layer can only exist where the measured sites share
genetic architecture in the first place. The six skeletal sites correlate at a mean of 0.44 with
each other; the seven brain volumes correlate at 0.12 (Figure 2C). Not one gene reaches Z > 3 at
all seven volumes, so there is no shared layer in the brain to detect. This is the phenome-wide
null result<sup>2</sup> seen from inside: averaged across organ systems, most of
which do not offer anatomically matched measurements on both sides, the correspondence we observe
in the skeleton is diluted until it disappears.

Everything so far is derived from one cohort, so we asked whether the same genes behave the same
way elsewhere. GEFOS is a meta-analysis of bone mineral density at the femoral neck, lumbar spine
and forearm that contains no UK Biobank participants.<sup>5</sup> Breadth defined in
UK Biobank reproduced the dose-response curve there without modification (Figure 3A): the median
gene-level Z in GEFOS rises monotonically with breadth, and shared-layer genes reach a median of
2.43 at the femoral neck and 2.12 at the lumbar spine against 0.11 and 0.08 for all other genes
(P = 4 × 10⁻³⁹ and 3 × 10⁻³²). Breadth survives the strictest available adjustment: entered
alongside all six UK Biobank site scores, from which it is derived, it retains t = 11.5 and 12.4.
Breadth is not a summary of the six scores; it is information they do not contain individually.

The same property reaches the clinic. In a fracture genome-wide association study of 426,795
individuals,<sup>6</sup> the proportion of genes with a fracture signal at Z > 2
rises from 5% at breadth zero to 53% at breadth six (Figure 3B), and breadth again holds beside
all six site scores (t = 11.2, P = 3 × 10⁻²⁹). Across replication and fracture alike, the
site-specific layer behaves like the background rather than like the shared layer: its median
GEFOS Z at the lumbar spine is 0.11 against 2.12 for the shared layer, and its fracture median is
0.97 against 2.33 (Figure 3C). Acting strongly at one skeletal site does not make a gene a
skeletal disease gene.

What is the shared layer made of? It is a coherent developmental module rather than a list of
statistical outliers (Figure 4A). Against all other genes, shared-layer genes are enriched for
FGF signaling (6 of 103, OR 16.8, P < 10⁻⁴), ossification (12 of 103, OR 6.4), canonical Wnt
signaling (9 of 103, OR 6.3), Wnt signaling overall (11 of 103, OR 5.4) and cartilage development
(5 of 103, OR 5.4), while BMP, Hedgehog and TGF-β show nothing and Notch falls short of
significance (3 of 103, OR 3.7, P = 0.07). Curated pathways carry the
literature's own biases, so we added a measure that is counted rather than curated: the density of
open chromatin around each gene's regulatory elements in human embryonic limb cartilage, at the
proximal and distal end of each of four long bones.<sup>7</sup> Shared-layer
genes sit in denser regulatory neighborhoods than length-matched expectation (13.2 against
11.1 ± 0.6 peaks per element, Z = +3.7, P = 2 × 10⁻⁴), and site-specific genes do not (Z = +0.2,
Figure 4B). The layer is the Wnt and ossification axis of skeletal development, and it is the
axis that contains the textbook of osteoporosis genetics — ESR1, LRP5, WNT16, WNT4, WLS, the
OPG–RANK pair, SOX6, CPED1, FGFRL1, IDUA, MDK and NOTUM (Figure 4C).

That the skeleton has site-specific genetics is not new. Kemp and colleagues showed a decade ago
that bone mineral density loci differ between skeletal sites, with the skull standing
apart,<sup>8</sup> and skull-specific signals have since been mapped in their own
right.<sup>9</sup> What the anatomy of bone mineral density has not been used
for is the disease side. Site-specificity has been read as a statement about which locus matters
where; here it is the complement, the set of genes that matter everywhere, that turns out to carry
the Mendelian burden and the fracture risk, while the site-specific genes carry neither.

One honest qualification governs how far this can be pushed. Within the shared layer, fracture
signal and Mendelian disease are not associated with each other: 33 of the 103 genes carry a
Mendelian skeletal condition, 55 carry a fracture signal, and 20 carry both, which is what
independence predicts. The defensible statement is that the two clinical ends draw on the same
layer, not that they draw on the same genes. A second qualification concerns transfer. When the
layer is redefined from scratch in GEFOS rather than imported from UK Biobank, the resulting gene
lists overlap by only 19% and the Mendelian enrichment weakens. The property replicates; the
particular list does not. We therefore report breadth throughout as a continuous property of a
gene, and treat the 103-gene list as an illustration rather than a deliverable.

The practical consequence is a prior that costs nothing to apply. A gene of unknown significance
that acts across the whole skeleton is a better candidate for a generalized skeletal dysplasia
than one that acts at a single site, and the prior can be read off a table that already exists.
More generally, the result tells the field where to look for the rare-common correspondence: not
in counts of traits across the phenome, but inside single organ systems where the same anatomy can
be named on both sides. The skeleton permits that naming. Where an organ does not — where its
measurement sites do not share genetic architecture, or where the ontology describes the organ in
a different vocabulary than the assay — the correspondence will be real and still invisible.

---

## Material and methods

**Genome-wide association data.** Bone mineral density at six DXA-derived skeletal sites was taken
from UK Biobank summary statistics released through the GWAS Catalog:<sup>3</sup> head
(GCST90568450, n = 31,986), arms (GCST90568449, n = 31,873), legs (GCST90568451, n = 31,873),
pelvis (GCST90568452, n = 31,873), femoral neck (GCST90568458, n = 32,017) and lumbar spine
(GCST90568448, n = 30,449). The same release contains five further DXA regions; the six above are
the ones with distinct anatomical names that the phenotype ontology also names. Heel bone mineral
density (GCST006979, n = 426,824) and fracture (GCST006980, 53,184 cases and 373,611 controls)
were taken from a published osteoporosis meta-analysis.<sup>6</sup> Replication used GEFOS/UK10K
femoral neck, lumbar spine and forearm,<sup>5</sup> which share no participants with UK Biobank;
those files carry no sample-size column, so the per-trait N required by the gene-level method was
estimated from standard errors and effect-allele frequencies as N = median[1 / (2f(1 − f)se²)]
over variants with 0.1 < f < 0.9, giving 34,806, 25,759 and 8,112. The cross-organ control used
seven subcortical volumes (thalamus, caudate, putamen, pallidum, hippocampus, amygdala,
accumbens) from Oxford BIG40, up to n = 33,211 each.<sup>4</sup>

**Gene-level statistics.** Variant-level summary statistics were harmonized to GRCh37, mapped to
genes with a 35 kb upstream and 10 kb downstream window, and converted to gene-level Z statistics
with MAGMA v1.10 using the 1000 Genomes European reference panel.<sup>10</sup> MAGMA
reports a p-value floor at 5 × 10⁻¹⁰; all analyses therefore use rank-based statistics rather than
means of the transformed values. Genes were retained if a Z statistic was available at all six
skeletal sites, leaving 18,392 genes.

**Breadth, the shared layer and the site-specific layer.** Breadth is the number of the six sites
at which a gene reaches Z > 2. The shared layer is the set of genes with breadth 6. The
site-specific layer is the set of genes reaching Z > 4 at one site while exceeding every other
site by more than 2 Z units, excluding any gene already in the shared layer. Threshold sensitivity
was assessed by repeating the definition at Z > 2.5, 3.0, 3.5 and 4.0.

**Locus-level analysis.** To remove linkage disequilibrium and paralogue clusters, genes were
sorted by position and merged into loci whenever consecutive genes on the same chromosome lay
within a given distance; merge distances of 0, 100, 250, 500 and 1,000 kb were used. A locus was
called shared-layer or Mendelian if any constituent gene was.

**Mendelian truth sides.** Three curations were used independently. Genomics England PanelApp
panel 309 (skeletal dysplasia) contributed green-rated genes.<sup>11</sup> The Human Phenotype
Ontology contributed genes annotated below the skeletal system branch, and, for the
phenotype-breadth analysis, the six body regions (skull, spine, thorax, pelvis, long bones, hands
and feet) each gene's terms touch.<sup>12</sup> ClinVar contributed genes with at least one
pathogenic or likely pathogenic variant whose condition field names a skeletal diagnosis, with
benign variants retained as a negative control.<sup>13</sup> The brain truth side was built by the
same procedure from the nervous-system branch of the ontology.

**Statistics.** Enrichments are Fisher exact tests with Haldane-corrected odds ratios and Wald
95% confidence intervals; unless stated otherwise the comparison group is every gene outside the
layer being tested. Dose-response trends are logistic regressions of the curation indicator
on breadth as a continuous variable. Phenotype breadth was analyzed by ordinary least squares on
breadth with the logarithm of the gene's total ontology term count as a covariate, and by
Mann-Whitney comparison of the extreme breadth strata. Conditional tests enter breadth alongside
all six site-level Z statistics in one linear model. Regulatory density was compared against a
null distribution of 5,000 resamples matched on the decile of log gene length. All tests are
two-sided unless a direction was specified in advance. The analysis was exploratory and not
preregistered; no correction for multiple comparisons is applied across the independent questions.
The threshold sweep, the locus sweep, the matched-null confounder controls and the full pathway
table are given in the supplemental information.

**Pathways and regulatory density.** Gene sets for FGF, Wnt, canonical Wnt, BMP, Hedgehog, TGF-β
and Notch signaling, ossification and cartilage development were taken from Gene Ontology
biological process terms, using non-electronic evidence codes only.<sup>14</sup> Open-chromatin
peak density was computed within 100 kb of each gene and expressed per regulatory element, from
ATAC-seq peaks called on hg19 in human embryonic limb cartilage at eight skeletal elements — the
proximal and distal end of the femur, tibia, humerus and radius, at embryonic day 54
(GEO: GSE252289).<sup>7</sup>

---

## Acknowledgments

To be completed.

## Author contributions

A.S. designed the study, wrote the analysis code, performed the analyses, produced the figures
and wrote the manuscript. C.T.T. supervised the work and revised the manuscript. Both authors
read and approved the final version.

## Declaration of interests

The authors declare no competing interests.

## Data and code availability

All data analyzed in this study are publicly available under the accessions listed in the
Material and methods: the summary statistics through the GWAS Catalog, Oxford BIG40 and
gefos.org, and the open-chromatin peaks through GEO accession GSE252289. The code that reproduces
every number and figure in this paper, together with the derived result tables from which each
figure is drawn, is available at https://github.com/alexschulzcell/skeletal-paper.

## Web resources

Human Phenotype Ontology, https://hpo.jax.org
PanelApp, https://panelapp.genomicsengland.co.uk
ClinVar, https://www.ncbi.nlm.nih.gov/clinvar
GWAS Catalog, https://www.ebi.ac.uk/gwas
Oxford BIG40, https://open.win.ox.ac.uk/ukbiobank/big40
MAGMA, https://cncr.nl/research/magma

## Declaration of generative AI and AI-assisted technologies in the manuscript preparation process

During the preparation of this work the authors used a generative AI assistant (Anthropic Claude)
to draft and review analysis and figure code, to typeset the figures, and to edit language. After
using this tool the authors reviewed and edited the content as needed and take full responsibility
for the content of the publication. No text, figure or result was accepted without verification
against the underlying data.

---

## References

1. Freund, M.K., Burch, K.S., Shi, H., Mancuso, N., Kichaev, G., Garske, K.M., Pan, D.Z., Miao,
   Z., Mohlke, K.L., Laakso, M., et al. (2018). Phenotype-specific enrichment of Mendelian
   disorder genes near GWAS regions across 62 complex traits. Am. J. Hum. Genet. *103*, 535–552.
2. Barbitoff, Y.A., Bogaichuk, P.M., Pavlova, N.S., Malysheva, P.V., and Predeus, A.V. (2025).
   Functional determinants and evolutionary consequences of pleiotropy in complex and Mendelian
   traits. Mol. Biol. Evol. *42*, msaf232.
3. Qian, Y., Xia, J., Wang, P., Xie, C., Lin, H.L., Li, G.H., Yuan, C.D., Qiu, M.C., Fang, Y.H.,
   Yu, C.F., et al. (2025). Genome-wide association studies of over 30,000 samples with bone
   mineral density at multiple skeletal sites and its clinical relevance. Genomics Proteomics
   Bioinformatics *23*, qzaf097.
4. Smith, S.M., Douaud, G., Chen, W., Hanayik, T., Alfaro-Almagro, F., Sharp, K., and Elliott,
   L.T. (2021). An expanded set of genome-wide association studies of brain imaging phenotypes in
   UK Biobank. Nat. Neurosci. *24*, 737–745.
5. Zheng, H.F., Forgetta, V., Hsu, Y.H., Estrada, K., Rosello-Diez, A., Leo, P.J., Dahia, C.L.,
   Park-Min, K.H., Tobias, J.H., Kooperberg, C., et al. (2015). Whole-genome sequencing identifies
   EN1 as a determinant of bone density and fracture. Nature *526*, 112–117.
6. Morris, J.A., Kemp, J.P., Youlten, S.E., Laurent, L., Logan, J.G., Chai, R.C., Vulpescu, N.A.,
   Forgetta, V., Kleinman, A., Mohanty, S.T., et al. (2019). An atlas of genetic influences on
   osteoporosis in humans and mice. Nat. Genet. *51*, 258–266.
7. Richard, D., Muthuirulan, P., Young, M., Yengo, L., Vedantam, S., Marouli, E., Bartell, E.,
   GIANT Consortium, Hirschhorn, J., and Capellini, T.D. (2025). Functional genomics of human
   skeletal development and the patterning of height heritability. Cell *188*, 15–32.e24.
8. Kemp, J.P., Medina-Gomez, C., Estrada, K., St Pourcain, B., Heppe, D.H.M., Warrington, N.M.,
   Oei, L., Ring, S.M., Kruithof, C.J., Timpson, N.J., et al. (2014). Phenotypic dissection of
   bone mineral density reveals skeletal site specificity and facilitates the identification of
   novel loci in the genetic regulation of bone mass attainment. PLoS Genet. *10*, e1004423.
9. Medina-Gomez, C., Mullin, B.H., Chesi, A., Prijatelj, V., Kemp, J.P., Shochat-Carvalho, C.,
   Trajanoska, K., Wang, C., Joro, R., Evans, T.E., et al. (2023). Bone mineral density loci
   specific to the skull portray potential pleiotropic effects on craniosynostosis. Commun. Biol.
   *6*, 691.
10. de Leeuw, C.A., Mooij, J.M., Heskes, T., and Posthuma, D. (2015). MAGMA: generalized gene-set
    analysis of GWAS data. PLoS Comput. Biol. *11*, e1004219.
11. Martin, A.R., Williams, E., Foulger, R.E., Leigh, S., Daugherty, L.C., Niblock, O., Leong,
    I.U.S., Smith, K.R., Gerasimenko, O., Haraldsdottir, E., et al. (2019). PanelApp crowdsources
    expert knowledge to establish consensus diagnostic gene panels. Nat. Genet. *51*, 1560–1565.
12. Gargano, M.A., Matentzoglu, N., Coleman, B., Addo-Lartey, E.B., Anagnostopoulos, A.V.,
    Anderton, J., Avillach, P., Bagley, A.M., Bakštein, E., Balhoff, J.P., et al. (2024). The
    Human Phenotype Ontology in 2024: phenotypes around the world. Nucleic Acids Res. *52*,
    D1333–D1346.
13. Landrum, M.J., Chitipiralla, S., Brown, G.R., Chen, C., Gu, B., Hart, J., Hoffman, D., Jang,
    W., Kaur, K., Liu, C., et al. (2020). ClinVar: improvements to accessing data. Nucleic Acids
    Res. *48*, D835–D844.
14. Gene Ontology Consortium, Aleksander, S.A., Balhoff, J., Carbon, S., Cherry, J.M., Drabkin,
    H.J., Ebert, D., Feuermann, M., Gaudet, P., Harris, N.L., et al. (2023). The Gene Ontology
    knowledgebase in 2023. Genetics *224*, iyad031.

---

## Figure titles and legends

**Figure 1. Breadth of action across the skeleton marks the Mendelian layer**

(A) Bone mineral density was measured at six anatomically named sites in UK Biobank. Each gene
receives one number, its breadth, defined as the count of sites at which the gene-level statistic
exceeds Z > 2. Filled circles illustrate a gene acting at all six sites and a gene acting at one.
(B) Proportion of genes carrying a Mendelian skeletal disease annotation as a function of breadth,
expressed relative to genes acting at no site, in three independent curations. Shaded bands are
95% Jeffreys intervals. Absolute proportions at breadth 6 and the per-site odds ratio from
logistic regression are given in the panel. n = 18,392 genes.
(C) Odds ratio for Mendelian skeletal disease when the shared layer is redefined at Z thresholds
from 2.0 to 4.0, which reduces it from 103 to 14 genes. The comparison group is every gene outside
the layer. Error bars are 95% confidence intervals.
(D) Odds ratio at the locus level, after collapsing neighboring genes into single loci at merge
distances from zero to one megabase. Error bars are 95% confidence intervals.

**Figure 2. The rule relates anatomy to anatomy and stops at the organ boundary**

(A) Mean number of the six body regions affected by the Mendelian phenotype, as a function of
breadth of action, for 2,789 genes carrying a skeletal condition. Error bars are standard errors
of the mean. The straight line is the fitted trend: the regression coefficient of 0.078 regions
per site, adjusted for the logarithm of each gene's total ontology term count, drawn through the
sample means.
(B) Cross-organ control. Each cell is the coefficient of breadth in an ordinary least squares
model of phenotype breadth, adjusted for ontology term count; cell shading is the t statistic.
Skeletal breadth predicts skeletal phenotype breadth; brain breadth predicts neither outcome.
(C) Mean genetic correlation between measurement sites within each organ. Bars are the mean across
all pairs; white lines span the minimum and maximum pairwise correlation. A shared layer requires
sites that share genetic architecture.

**Figure 3. The layer replicates in an independent cohort and carries fracture risk**

(A) Median gene-level Z in GEFOS, a cohort with no UK Biobank participants, as a function of
breadth defined in UK Biobank, for three skeletal sites. The quoted P value is the shared layer
against all other genes at the femoral neck.
(B) Percentage of genes with a fracture-risk signal at Z > 2 as a function of breadth. The
conditional test enters breadth alongside all six UK Biobank site scores.
(C) Median gene-level Z for the shared layer, the site-specific layer and all other genes, across
two replication endpoints and fracture. n = 103 and 36 genes in the two layers, against 18,245
other genes for the GEFOS endpoints and 18,247 for fracture.

**Figure 4. What the shared layer is made of**

(A) Enrichment of the shared layer for developmental signaling pathways, against all other genes.
Filled circles mark P < 0.05; bars are 95% confidence intervals; counts give the number of
shared-layer genes in each pathway.
(B) Density of open chromatin in human fetal limb tissue around the regulatory elements of each
layer. The grey bar and whiskers are the mean and 95% interval of 5,000 length-matched resamples;
the filled circle is the observed value.
(C) The 103 shared-layer genes, positioned by their fracture-risk signal and by their independent
replication in GEFOS. Filled symbols carry a Mendelian skeletal disease annotation. Dashed lines
mark Z = 2 on each axis. The labeled genes are those named in the main text.

**Graphical abstract.** A gene acting at every skeletal site produces generalized skeletal
dysplasia when severely mutated and low bone density with fracture when finely tuned by common
variation. A gene acting at a single site does neither.
