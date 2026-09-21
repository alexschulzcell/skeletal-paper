#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Gene-level association statistics for every trait, with MAGMA.
#
# What it does   annotates variants to genes once (35 kb upstream, 10 kb
#                downstream), then runs MAGMA's gene analysis for each trait in
#                turn, several at a time. Every trait is independent, so the
#                script is restartable and parallelisable: a trait whose
#                `.genes.out` already exists is skipped.
# Reads          $SKELBREADTH_DATA/magma/<trait>.pval        (stage 01, script 01)
#                $SKELBREADTH_DATA/magma/sample_sizes.tsv    (stage 01, script 01)
#                $SKELBREADTH_DATA/reference/g1000_eur.*
#                $SKELBREADTH_DATA/reference/NCBI37.3.gene.loc
# Writes         $SKELBREADTH_DATA/magma/annot.genes.annot
#                $SKELBREADTH_DATA/magma/<trait>.genes.out
#                $SKELBREADTH_DATA/magma/<trait>.genes.raw
#                results/magma_run.log
# Numbers        none directly. Every Z statistic in the paper descends from
#                these files.
#
# Runtime        THIS IS THE LONG STEP. Two to six hours for the full set of
#                21 traits on a four-core laptop, longer for the height file.
#                It is not hung. Progress is written to results/magma_run.log.
# Hardware       8 GB RAM is enough; 16 GB if you raise -j above 4. About
#                20 GB of disk for the .genes.raw files, which stage 02 does
#                not need but which are kept so that a gene-set analysis can
#                be added later without recomputing.
#
# The annotation window is set in 00_setup/config.py (MAGMA_WINDOW). The paper
# reports 35/10 kb and a sensitivity analysis at 10/10 kb; see
# docs/KNOWN_ISSUES.md. To produce the sensitivity set, run:
#     MAGMA_WINDOW=10,10 bash 01_gene_scores/02_run_magma.sh --suffix _w10
# ---------------------------------------------------------------------------
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${SKELBREADTH_ROOT:-$(dirname "$HERE")}"
DATA="${SKELBREADTH_DATA:-$ROOT/data_raw}"
MG="$DATA/magma"
REF="$DATA/reference"
MAGMA="${MAGMA_BIN:-magma}"
WINDOW="${MAGMA_WINDOW:-35,10}"
JOBS="${JOBS:-4}"
SUFFIX=""

while [ $# -gt 0 ]; do
  case "$1" in
    --suffix) SUFFIX="$2"; shift 2 ;;
    -j) JOBS="$2"; shift 2 ;;
    *) echo "unknown argument: $1"; exit 2 ;;
  esac
done

mkdir -p "$ROOT/results"
LOG="$ROOT/results/magma_run.log"

log() { echo "$(date '+%Y-%m-%d %H:%M:%S')  $*" | tee -a "$LOG"; }

# --- preflight -------------------------------------------------------------
command -v "$MAGMA" >/dev/null 2>&1 || {
  echo "MAGMA not found. Set MAGMA_BIN, or install it:"
  echo "  conda env create -f environment.yml && conda activate skeletal-breadth"
  exit 1; }
[ -s "$REF/g1000_eur.bim" ] || { echo "MISSING: $REF/g1000_eur.bim - run 00_setup/02_download_references.sh"; exit 1; }
[ -s "$REF/NCBI37.3.gene.loc" ] || { echo "MISSING: $REF/NCBI37.3.gene.loc"; exit 1; }
[ -s "$MG/sample_sizes.tsv" ] || { echo "MISSING: $MG/sample_sizes.tsv - run 01_gene_scores/01_harmonise_sumstats.py"; exit 1; }

log "MAGMA:  $("$MAGMA" --version 2>&1 | head -1)"
log "window: $WINDOW kb   jobs: $JOBS"

# --- annotation, once for all traits ---------------------------------------
ANNOT="$MG/annot${SUFFIX}"
if [ ! -s "$ANNOT.genes.annot" ]; then
  log "annotating variants to genes (window $WINDOW kb) ..."
  "$MAGMA" --annotate window="$WINDOW" \
           --snp-loc "$REF/g1000_eur.bim" \
           --gene-loc "$REF/NCBI37.3.gene.loc" \
           --out "$ANNOT" >> "$LOG" 2>&1 \
    || { log "annotation FAILED - see $LOG"; exit 1; }
  log "annotation done: $(wc -l < "$ANNOT.genes.annot") genes"
else
  log "annotation already present, reused"
fi

# --- one gene analysis per trait -------------------------------------------
# NOTE: MAGMA appends its own extensions to --out. Do not pass a filename that
# already ends in .genes.out; we chased a non-existent "path bug" for an hour
# over exactly that.
run_trait() {
  local trait="$1" n="$2"
  if [ -s "$MG/${trait}${SUFFIX}.genes.out" ]; then
    log "  $trait: already done, skipped"
    return 0
  fi
  [ -s "$MG/$trait.pval" ] || { log "  $trait: no .pval file, skipped"; return 0; }
  log "  $trait: starting (N $n)"
  "$MAGMA" --bfile "$REF/g1000_eur" \
           --pval "$MG/$trait.pval" N="$n" \
           --gene-annot "$ANNOT.genes.annot" \
           --out "$MG/${trait}${SUFFIX}" > "$MG/${trait}${SUFFIX}.run.log" 2>&1
  if [ -s "$MG/${trait}${SUFFIX}.genes.out" ]; then
    log "  $trait: done ($(( $(wc -l < "$MG/${trait}${SUFFIX}.genes.out") - 1 )) genes)"
  else
    log "  $trait: FAILED - see $MG/${trait}${SUFFIX}.run.log"
  fi
}

log "=== gene analysis ==="
running=0
# tail -n +2 skips the header; the sample-size file is written newline-clean by
# the harmonisation step, so no \r can reach MAGMA's N= argument here.
while IFS=$'\t' read -r trait n source; do
  [ -z "${trait:-}" ] && continue
  run_trait "$trait" "$n" &
  running=$(( running + 1 ))
  if [ "$running" -ge "$JOBS" ]; then wait -n 2>/dev/null || wait; running=$(( running - 1 )); fi
done < <(tail -n +2 "$MG/sample_sizes.tsv")
wait

log "=== summary ==="
done_n=0; fail=()
while IFS=$'\t' read -r trait n source; do
  [ -z "${trait:-}" ] && continue
  if [ -s "$MG/${trait}${SUFFIX}.genes.out" ]; then done_n=$(( done_n + 1 ));
  else fail+=("$trait"); fi
done < <(tail -n +2 "$MG/sample_sizes.tsv")
log "$done_n trait(s) complete"
if [ ${#fail[@]} -gt 0 ]; then
  log "FAILED: ${fail[*]}"
  log "The six DXA sites are required; the endpoints are not."
fi

echo
echo "Now run:  python 01_gene_scores/03_assemble_gene_matrix.py"
