# Einreichung — AJHG, Format „Report"

Stand 21.09.2026. Alle Vorgaben aus cell.com/ajhg/authors, abgerufen am 21.09.2026.

Diese Datei gehört **nicht** in den Einreichungsordner. Sie bleibt hier, im
Repository. `06_manuscript/assemble_submission.py` leert `../submission/` und
legt dort ausschließlich die acht Dateien ab, die die Zeitschrift bekommt.

## Was eingereicht wird

| Datei | Zweck | Quelle |
|---|---|---|
| `Manuscript.pdf` | Haupttext, Zeilen- und Seitennummern, Abbildungen eingebettet | `manuscript.md` → `render_pdf.py` |
| `SupplementalInformation.pdf` | 13 Tabellen und zwei Anmerkungen | `build_supplement.py` → `render_pdf.py` |
| `CoverLetter.pdf` | vertraulich, Reviewer sehen ihn nicht | `cover_letter.md` → `render_pdf.py` |
| `Figure1–4.tif` | 174 mm breit, 600 dpi, LZW, RGB ohne Alphakanal | `05_figures/fig*.py` |
| `GraphicalAbstract.tif` | 1650 × 1650 px = 5,5 Zoll bei 300 dpi | `05_figures/graphical_abstract.py` |

Nichts davon wird von Hand bearbeitet. Der ganze Ordner entsteht aus

```bash
make figures && make manuscript
python 06_manuscript/assemble_submission.py
```

Die Vektorfassungen der Abbildungen (`figures/*.pdf`) bleiben im Repository:
sie sind für die eigene Durchsicht, nicht für die Einreichung.

## Vorgaben und wie wir sie erfüllen

| Vorgabe | Wert | erfüllt |
|---|---|---|
| Report, ca. 3.600 Wörter | 2.589 (Haupttext 1.895 + Methoden 694) | ✓ |
| Abstract ≤ 250 Wörter, ein Absatz, ohne Zitate | 250 | ✓ |
| Titel ≤ 3 Zeilen à 54 Zeichen | 2 Zeilen, 42 und 47 Zeichen | ✓ |
| Haupttext ohne Untergliederung | durchlaufend, Methoden separat (bei Reports erlaubt) | ✓ |
| Abbildungsbreite 85 / 114 / 174 mm | alle exakt 174,00 mm, von `style.save()` erzwungen | ✓ |
| Schrift in Abbildungen ≥ 6 pt | Arial, 6–9 pt, ohne nachträgliche Skalierung | ✓ |
| Eine Schriftfamilie pro Abbildung | nur Arial, auch in den Exponenten | ✓ |
| TIFF ohne Alphakanal, LZW, ≥ 600 dpi | RGB, LZW, exakt 600 dpi | ✓ |
| Panel-Buchstaben Großbuchstaben, fett, oben links | ✓ | ✓ |
| Legende: Titel ohne Panel-Bezug, dann Panel für Panel | ✓ | ✓ |
| Graphical Abstract: ein Panel, 5,5 Zoll, 300 dpi, Arial 12–16 pt | ✓ | ✓ |
| GA ohne Datenelemente, liest oben → unten | ✓ | ✓ |
| GA nicht mit Allzweck-Bildgenerator erzeugt | vollständig aus matplotlib-Code | ✓ |
| Zitate als hochgestellte Ziffern, in Reihenfolge | 14 Referenzen, geprüft gegen PubMed | ✓ |
| „et al." erst ab elf Autoren | erste zehn, dann „et al." | ✓ |
| Declaration of interests | vorhanden | ✓ |
| KI-Erklärung vor der Referenzliste | vorhanden | ✓ |
| Zeilen- und Seitennummern für die Begutachtung | durchlaufend, im PDF | ✓ |
| Erstfassung als ein PDF < 15 MB | 1,7 MB, Abbildungen eingebettet | ✓ |
| Data and code availability | vollständig, Repo öffentlich | ✓ |
| Autorenliste, Affiliation, ORCIDs, Correspondence | vollständig | ✓ |
| Acknowledgments | offen | ☐ |

## Was vor dem Absenden noch zu tun ist

Drei Punkte, die niemand außer den Autoren schließen kann.

1. **Abbildung 4B neu rechnen.** Die ATAC-Akzession im Repository war falsch
   (GSE170199 = HepG2-ChIP-seq). Die richtige Quelle ist identifiziert und
   verdrahtet — GSE252289, Richard et al., *Cell* 2025, acht Limb-Elemente,
   hg19 — aber die Zahlen in Abbildung 4B stammen noch aus dem alten Lauf mit
   Peak-Dateien unbekannter Herkunft.

   Eine Teilkontrolle gegen GSE252289 läuft bereits ohne Pipeline:

   ```bash
   python 00_setup/04_prepare_atac_peaks.py
   python 04_analysis/07b_recheck_regulatory_density.py
   ```

   Ergebnis: 7,50 gegen einen längengematchten Null von 5,50 ± 0,46,
   Z = +4,39 — gegenüber 13,16 / 11,11 ± 0,55 / Z = +3,69 in der
   ausgelieferten Tabelle. **Die Aussage hält, die absoluten Werte nicht.**
   Der Überschuss über den Null ist relativ sogar größer (1,36-fach statt
   1,18-fach).

   Die Teilkontrolle kann das Panel nicht ersetzen: sie nutzt als Hintergrund
   alle Gene der MAGMA-Positionstabelle statt der 18.392 mit Z an allen sechs
   Sites, und die 36 site-spezifischen Gene — der zweite Punkt im Panel —
   lassen sich aus dem Ausgelieferten gar nicht rekonstruieren. Für die
   endgültige Fassung:

   ```bash
   python 00_setup/04_prepare_atac_peaks.py
   python 04_analysis/07_pathways_and_axis.py    # braucht die volle Pipeline
   make figures && make verify
   ```

   Danach die drei Zahlen im Ergebnisteil („13.2 against 11.1 ± 0.6 …,
   Z = +3.7, P = 2 × 10⁻⁴") und in `numbers.json` anpassen. Siehe
   `docs/KNOWN_ISSUES.md`, Punkt 12.
2. **Acknowledgments** ausfüllen; steht als „To be completed".
3. **Reviewer-Vorschläge und Ausschlüsse** in den Cover Letter.

## Kontrolle vor jedem Export

```bash
make verify
```

Prüft 58 Zahlen im Manuskript gegen die Tabellen in `repo/results/` und bricht
bei jeder Abweichung ab. Muss grün sein, bevor etwas das Haus verlässt.

Die Abbildungsskripte prüfen sich ebenfalls selbst: `05_figures/style.py`
weigert sich, eine Abbildung zu schreiben, deren Leinwand nicht exakt die
verlangte Spaltenbreite hat, deren Inhalt über den Rand läuft oder in der sich
zwei Textstücke überlappen.

## Zwei Zahlen ohne Tabelle

Zwei Angaben im Haupttext stehen in keiner Datei unter `results/`: die
Korrelation von −0,013 zwischen Breite und Variantenzahl pro Gen, und die
Aussage, dass benigne ClinVar-Varianten keine Anreicherung zeigen. Beide
stammen aus `04_analysis/04_confounders.py` bzw. `01_enrichment.py`, deren
Ausgaben nicht zu den Panel-Tabellen gehören und deshalb von `make verify`
nicht erfasst werden. Vor der Einreichung einmal gegen
`results/confounder_matched_nulls.tsv` und `results/enrichment_layers.tsv`
gegenprüfen.
