#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Download every GWAS summary statistic the paper uses.
#
# What it does   fetches the discovery traits (six DXA sites, heel BMD,
#                fracture, otosclerosis, height, sitting-height ratio), the
#                seven subcortical brain volumes used for the cross-organ
#                control, and the three GEFOS replication traits.
# Reads          nothing; all sources are public archives.
# Writes         $SKELBREADTH_DATA/sumstats/<trait>.tsv.gz
#                $SKELBREADTH_DATA/sumstats/md5sums.txt
# Numbers        none. This script produces no result; it produces the inputs
#                of 01_gene_scores/.
#
# Size and time  roughly 40 GB and, on a 100 Mbit line, three to six hours.
#                The script is restartable: a trait whose file is already
#                present and non-empty is skipped.
#
# One trait cannot be downloaded: otosclerosis (N 864,702) is our own
# harmonisation of an in-house case/control analysis and is not in a public
# archive. See docs/DATA_SOURCES.md. Every otosclerosis result in the paper is
# marked as such, and 04_analysis/06_clinical_endpoints.py skips that endpoint
# with a stated message when the file is absent rather than failing.
# ---------------------------------------------------------------------------
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${SKELBREADTH_ROOT:-$(dirname "$HERE")}"
DATA="${SKELBREADTH_DATA:-$ROOT/data_raw}"
OUT="$DATA/sumstats"
mkdir -p "$OUT"

GWAS_FTP="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
BIG40="https://open.win.ox.ac.uk/ukbiobank/big40/release2/stats33k"

# trait  accession/id  url
# The six DXA sites are the Human Genetics 2023 UK Biobank DXA release,
# GCST90568448-GCST90568458; the three we do not use are skipped.
read -r -d '' TRAITS <<'EOF' || true
head     GCST90568448  ebi
arms     GCST90568449  ebi
legs     GCST90568450  ebi
pelvis   GCST90568451  ebi
fn       GCST90568452  ebi
ls       GCST90568453  ebi
heel     GCST90568454  ebi
frak     GCST006980    ebi
prop     GCST90728588  ebi
EOF

fetch() {   # fetch <url> <destination>
  local url="$1" dest="$2"
  if [ -s "$dest" ]; then
    echo "  already present: $(basename "$dest")"
    return 0
  fi
  echo "  downloading $(basename "$dest")"
  curl -fL --retry 5 --retry-delay 10 -C - -o "$dest.part" "$url"
  mv "$dest.part" "$dest"
}

ebi_url() {  # EBI lays summary statistics out in buckets of 1000 accessions
  local acc="$1" num bucket_lo bucket_hi
  num="${acc#GCST}"
  bucket_lo=$(( (10#$num / 1000) * 1000 + 1 ))
  bucket_hi=$(( (10#$num / 1000 + 1) * 1000 ))
  printf '%s/GCST%08d-GCST%08d/%s/%s.tsv.gz\n' \
    "$GWAS_FTP" "$bucket_lo" "$bucket_hi" "$acc" "$acc"
}

echo "=== GWAS Catalog traits ==="
while read -r trait acc src; do
  [ -z "${trait:-}" ] && continue
  [ "$src" = "ebi" ] || continue
  fetch "$(ebi_url "$acc")" "$OUT/$trait.tsv.gz"
done <<< "$TRAITS"

echo "=== Standing height (Yengo 2022, GIANT) ==="
fetch "https://portals.broadinstitute.org/collaboration/giant/images/4/4e/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR.gz" \
      "$OUT/height.tsv.gz"

echo "=== Seven subcortical brain volumes (Oxford BIG40, N 33,211 each) ==="
# BIG40 IDP numbers for the FIRST-segmented subcortical volumes.
declare -A BRAIN=(
  [thalamus]=0095 [caudate]=0096 [putamen]=0097 [pallidum]=0098
  [hippocampus]=0099 [amygdala]=0100 [accumbens]=0101
)
for region in "${!BRAIN[@]}"; do
  fetch "$BIG40/${BRAIN[$region]}.txt.gz" "$OUT/$region.tsv.gz"
done

echo "=== GEFOS / UK10K replication cohort (Zheng 2015, no UK Biobank) ==="
GEFOS="http://www.gefos.org/sites/default/files"
fetch "$GEFOS/fn2stu.MAF0_.005.pos_.out_.gz" "$OUT/fn2stu.tsv.gz"
fetch "$GEFOS/ls2stu.MAF0_.005.pos_.out_.gz" "$OUT/ls2stu.tsv.gz"
fetch "$GEFOS/fa2stu.MAF0_.005.pos_.out_.gz" "$OUT/fa2stu.tsv.gz"

echo "=== checksums ==="
( cd "$OUT" && md5sum ./*.tsv.gz > md5sums.txt )
echo "  wrote $OUT/md5sums.txt"
echo
echo "Now run:  python 00_setup/03_check_inputs.py"
