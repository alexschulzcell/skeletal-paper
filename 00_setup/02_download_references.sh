#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Download the reference files: LD panel, gene annotation, and the three
# independent curations that form the truth side of the paper.
#
# What it does   fetches the 1000 Genomes European LD reference and the
#                NCBI37.3 gene location file that MAGMA needs, the Human
#                Phenotype Ontology and its gene/disease annotations, the
#                ClinVar variant summary, the Genomics England PanelApp panel
#                PA309 (skeletal dysplasia, confidence 3), gnomAD constraint,
#                the GO term memberships used for the pathway analysis, and the
#                fetal limb ATAC peaks used for the regulatory measure.
# Reads          nothing.
# Writes         $SKELBREADTH_DATA/reference/*
# Numbers        none.
#
# Size and time  roughly 4 GB, twenty minutes on a fast line. Restartable.
# ---------------------------------------------------------------------------
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${SKELBREADTH_ROOT:-$(dirname "$HERE")}"
DATA="${SKELBREADTH_DATA:-$ROOT/data_raw}"
REF="$DATA/reference"
mkdir -p "$REF"

fetch() {
  local url="$1" dest="$2"
  if [ -s "$dest" ]; then echo "  already present: $(basename "$dest")"; return 0; fi
  echo "  downloading $(basename "$dest")"
  curl -fL --retry 5 --retry-delay 10 -C - -o "$dest.part" "$url"
  mv "$dest.part" "$dest"
}

echo "=== MAGMA auxiliary files (CNCR, Vrije Universiteit Amsterdam) ==="
MAGMA_AUX="https://vu.data.surfsara.nl/index.php/s"
fetch "$MAGMA_AUX/VnME7VUOoW1dZ1L/download" "$REF/g1000_eur.zip"
fetch "$MAGMA_AUX/Pj2orwuF3ZeqwDo/download" "$REF/NCBI37.3.zip"
( cd "$REF" && unzip -n -q g1000_eur.zip && unzip -n -q NCBI37.3.zip )
# leaves g1000_eur.{bed,bim,fam,synonyms} and NCBI37.3.gene.loc in place

echo "=== Human Phenotype Ontology ==="
fetch "https://purl.obolibrary.org/obo/hp.obo"                          "$REF/hp.obo"
fetch "https://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa"          "$REF/phenotype.hpoa"
fetch "https://purl.obolibrary.org/obo/hp/hpoa/genes_to_disease.txt"    "$REF/genes_to_disease.txt"

echo "=== ClinVar (weekly VCF, GRCh37) ==="
fetch "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz" \
      "$REF/clinvar.vcf.gz"

echo "=== Genomics England PanelApp, panel 309 (skeletal dysplasia) ==="
if [ ! -s "$REF/panelapp_pa309.json" ]; then
  echo "  querying the PanelApp API"
  curl -fL --retry 5 \
    "https://panelapp.genomicsengland.co.uk/api/v1/panels/309/?format=json" \
    -o "$REF/panelapp_pa309.json"
else
  echo "  already present: panelapp_pa309.json"
fi

echo "=== gnomAD v2.1.1 constraint (LOEUF) ==="
fetch "https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz" \
      "$REF/gnomad_constraint.txt.bgz"

echo "=== Gene Ontology annotations (human, GOA) ==="
fetch "http://geneontology.org/gene-associations/goa_human.gaf.gz" "$REF/goa_human.gaf.gz"
fetch "https://purl.obolibrary.org/obo/go/go-basic.obo"            "$REF/go-basic.obo"

echo "=== Publication counts per gene (gene2pubmed, for the ascertainment control) ==="
fetch "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz" "$REF/gene2pubmed.gz"

echo "=== Human embryonic limb ATAC peaks (hg19), eight skeletal elements ==="
# GSE252289, the ATAC-seq series of Richard et al., Cell 2025. The peaks are
# published as one Excel workbook with a sheet per skeletal element, so the
# conversion to per-element BED needs Python and lives in its own script.
#
# This block used to fetch GSE170199 and call it fetal chromatin
# accessibility. GSE170199 is a two-sample HepG2 ChIP-seq series from ENCODE:
# wrong assay, wrong tissue. See docs/KNOWN_ISSUES.md, issue 12.
if ls "$REF"/peaks/*hg19.bed.gz >/dev/null 2>&1; then
  echo "  already present: $(ls "$REF"/peaks/*hg19.bed.gz | wc -l) peak files"
else
  "${PYTHON:-python}" "$(dirname "$0")/04_prepare_atac_peaks.py"
fi

echo
echo "Now run:  python 00_setup/03_check_inputs.py"
