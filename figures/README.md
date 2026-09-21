# `figures/`

The four main figures of the paper and the graphical abstract, as produced by
the scripts in `05_figures/`, in PDF (vector) and TIFF (submission raster).

Nothing in this directory is computed. Every figure script reads panel tables
from `results/` and draws them; if a number in a figure disagrees with the
corresponding table, the figure is wrong, not the table.

| file | panel content | built from |
|---|---|---|
| `figure1.*` | A the six sites and the definition of breadth; B the three-curation dose-response; C the threshold sweep; D the locus sweep | `fig1b_*.csv`, `fig1c_*.csv`, `fig1d_*.csv` |
| `figure2.*` | A phenotype breadth against breadth of action, with the ascertainment control; B the 2x2 cross-organ matrix; C inter-site genetic correlation | `fig2a_*.csv`, `fig2b_*.csv`, `fig2c_*.csv` |
| `figure3.*` | A the dose-response redrawn in GEFOS; B fracture risk by breadth; C the two layers across three endpoints | `fig3a_*.csv`, `fig3b_*.csv`, `fig3c_*.csv` |
| `figure4.*` | A pathway enrichment of the shared layer; B regulatory density against a length-matched null; C the 103 shared-layer genes across the two clinical ends | `fig4a_*.csv`, `fig4b_*.csv`, `fig4c_*.csv` |
| `graphical_abstract.*` | one schematic panel; no data | nothing |

## What the files are guaranteed to be

`05_figures/style.py` refuses to write a figure unless all of this holds, so
these are properties of the files and not intentions:

- the canvas is **exactly 174 mm** wide (5.5 inches square for the graphical
  abstract) and is never trimmed, so the type keeps the point size the script
  asked for;
- no ink falls outside the canvas;
- no two pieces of text overlap;
- the TIFF is **RGB without an alpha channel**, LZW-compressed, at 600 dpi
  (300 dpi for the graphical abstract);
- only one typeface is embedded, including in the mathematical exponents.

Figures are licensed CC BY 4.0; see `LICENSE-DATA`.
