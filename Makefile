# ===========================================================================
# skeletal-breadth - companion repository
#
# `make` on its own prints this help. Nothing runs until you ask for it,
# because two of the targets take hours.
#
# Stages run strictly in the order their directories are numbered. Each script
# reads only files that an earlier-numbered script wrote, so any stage can be
# re-run on its own without re-running the ones before it.
# ===========================================================================

PYTHON ?= python
BASH   ?= bash
JOBS   ?= 4

# Where the downloads land. Override to put them on another disk:
#   make data SKELBREADTH_DATA=/scratch/skelbreadth
export SKELBREADTH_DATA ?= $(CURDIR)/data_raw
export SKELBREADTH_ROOT ?= $(CURDIR)

.DEFAULT_GOAL := help
.PHONY: help data genes layers truth analysis figures manuscript verify all \
        clean clean-results check

# ---------------------------------------------------------------------------

help:
	@echo ""
	@echo "  skeletal-breadth - make targets"
	@echo "  ------------------------------------------------------------------"
	@echo "  make data       download all summary statistics and reference files,"
	@echo "                  then verify them            ~40 GB, 3-6 h"
	@echo "  make genes      harmonise, run MAGMA, assemble the gene matrix"
	@echo "                  THE LONG STEP               2-6 h, 8 GB RAM"
	@echo "  make layers     define breadth, loci and the threshold sweep    <1 min"
	@echo "  make truth      build the four truth sides                      ~10 min"
	@echo "  make analysis   layers + truth + the seven analyses             ~15 min"
	@echo "  make figures    redraw the figures from results/                ~2 min"
	@echo "  make manuscript supplement, then the three submission PDFs      ~1 min"
	@echo "  make all        data -> genes -> analysis -> figures -> paper"
	@echo ""
	@echo "  make check      verify the inputs without running anything      <1 min"
	@echo "  make verify     check every manuscript number against results/  <1 min"
	@echo "  make clean-results   delete results/ and figures/ contents"
	@echo "  make clean      also delete MAGMA intermediates (keeps downloads)"
	@echo ""
	@echo "  Variables:  JOBS=$(JOBS)   PYTHON=$(PYTHON)"
	@echo "              SKELBREADTH_DATA=$(SKELBREADTH_DATA)"
	@echo ""
	@echo "  From nothing to the finished figures is about eight hours, almost"
	@echo "  all of it in 'data' and 'genes'. Both are restartable: re-running"
	@echo "  them skips whatever is already complete."
	@echo ""

# --- stage 00 --------------------------------------------------------------

data:
	$(BASH) 00_setup/01_download_sumstats.sh
	$(BASH) 00_setup/02_download_references.sh
	$(PYTHON) 00_setup/03_check_inputs.py

check:
	$(PYTHON) 00_setup/03_check_inputs.py --quick

# --- stage 01 --------------------------------------------------------------

genes:
	$(PYTHON) 01_gene_scores/01_harmonise_sumstats.py
	JOBS=$(JOBS) $(BASH) 01_gene_scores/02_run_magma.sh
	$(PYTHON) 01_gene_scores/03_assemble_gene_matrix.py

# --- stage 02 --------------------------------------------------------------

layers:
	$(PYTHON) 02_layers/01_define_breadth.py
	$(PYTHON) 02_layers/02_locus_clumping.py
	$(PYTHON) 02_layers/03_threshold_sweep.py

# --- stage 03 --------------------------------------------------------------

truth:
	$(PYTHON) 03_truth_sides/01_build_hpo_regions.py
	$(PYTHON) 03_truth_sides/02_build_panelapp.py
	$(PYTHON) 03_truth_sides/03_build_clinvar.py
	$(PYTHON) 03_truth_sides/04_build_brain_truth.py

# --- stage 04 --------------------------------------------------------------

analysis: layers truth
	$(PYTHON) 04_analysis/01_enrichment.py
	$(PYTHON) 04_analysis/02_symmetry_law.py
	$(PYTHON) 04_analysis/03_cross_organ_control.py
	$(PYTHON) 04_analysis/04_confounders.py
	$(PYTHON) 04_analysis/05_replication_gefos.py
	$(PYTHON) 04_analysis/06_clinical_endpoints.py
	$(PYTHON) 04_analysis/07_pathways_and_axis.py

# --- stage 05 --------------------------------------------------------------
# Each script refuses to write a figure that is the wrong width, that has ink
# outside the canvas, or that contains overlapping text, so a failure here is
# a layout bug and not a warning.

figures:
	@for f in 05_figures/fig*.py 05_figures/graphical_abstract.py; do \
	  echo "--- $$f"; $(PYTHON) "$$f" || exit 1; \
	done

# --- stage 06 --------------------------------------------------------------

verify:
	$(PYTHON) 06_manuscript/check_numbers.py

manuscript: verify
	$(PYTHON) 06_manuscript/build_supplement.py
	$(PYTHON) 06_manuscript/render_pdf.py

all: data genes analysis figures manuscript
	@echo ""
	@echo "Done. Every number of the paper is in results/, every figure in"
	@echo "figures/, and the three submission PDFs in 06_manuscript/pdf/."

# --- housekeeping ----------------------------------------------------------

clean-results:
	@find results -type f ! -name '.gitkeep' ! -name 'README.md' \
	     ! -name 'fig*.csv' -delete
	@find figures -type f ! -name '.gitkeep' ! -name 'README.md' -delete
	@rm -rf 06_manuscript/pdf
	@echo "results/ and figures/ emptied, PDFs removed. Downloads untouched."
	@echo "The fig*.csv panel tables are kept: nothing in this repository"
	@echo "regenerates them yet. See issue 11 in docs/KNOWN_ISSUES.md."

clean: clean-results
	@rm -f $(SKELBREADTH_DATA)/magma/*.genes.raw \
	       $(SKELBREADTH_DATA)/magma/*.pval \
	       $(SKELBREADTH_DATA)/magma/*.log \
	       $(SKELBREADTH_DATA)/magma/*.run.log
	@echo "MAGMA intermediates removed. Downloads untouched; 'make genes' will"
	@echo "re-harmonise and re-run, which takes hours."
