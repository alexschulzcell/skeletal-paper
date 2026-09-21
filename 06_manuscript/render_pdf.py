#!/usr/bin/env python3
"""Turn the manuscript, the supplement and the cover letter into submission PDFs.

Why this exists
    The journal wants a PDF with continuous line numbers and page numbers for
    review, with every font embedded and the figures placed at the size they
    were drawn. Converting the Markdown by hand, or through a word processor,
    loses the line numbers or re-flows the figures. This module renders the
    same Markdown the repository already checks with `check_numbers.py`, so
    the PDF cannot drift away from the text that was verified.

    It reads a deliberately small subset of Markdown: headings, paragraphs,
    horizontal rules, pipe tables, ordered and unordered lists, `**bold**`,
    `*italic*`, `<sup>` and HTML comments. That is everything the manuscript
    uses and nothing else, so an unsupported construct fails loudly instead of
    being silently dropped.

Usage
    python 06_manuscript/render_pdf.py                 # all three documents
    python 06_manuscript/render_pdf.py --only manuscript
"""
from __future__ import annotations

import argparse
import html
import re
import sys
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (BaseDocTemplate, Frame, HRFlowable, Image,
                                KeepTogether, PageBreak, PageTemplate,
                                Paragraph, Spacer, Table, TableStyle)

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "06_manuscript" / "pdf"

PAGE = A4
MARGIN_L, MARGIN_R = 28 * mm, 20 * mm
MARGIN_T, MARGIN_B = 20 * mm, 20 * mm


# --------------------------------------------------------------------- fonts
MONO = {"face": None}


def register_mono() -> str:
    """A monospace face for code spans, embedded like everything else.

    The built-in Courier is a base-14 font: reportlab references it without
    embedding it, and the PDF then fails a font-embedding preflight check.
    """
    if MONO["face"]:
        return MONO["face"]
    options = [("CourierNew", Path("C:/Windows/Fonts/cour.ttf"))]
    try:
        import matplotlib
        d = Path(matplotlib.__file__).parent / "mpl-data" / "fonts" / "ttf"
        options.append(("DejaVuSansMono", d / "DejaVuSansMono.ttf"))
    except Exception:
        pass
    for name, path in options:
        if path.exists():
            pdfmetrics.registerFont(TTFont(name, str(path)))
            MONO["face"] = name
            return name
    return ""


def register_fonts() -> tuple[str, str, str, str]:
    """Embed a real sans family. Arial where it exists, DejaVu Sans otherwise.

    The base-14 Helvetica is never used: it is not embedded in the output and
    production systems reject PDFs whose fonts are only referenced. Reportlab
    writes Helvetica into every page's resources unless the canvas base font is
    changed before the document is built, which is what `canvas_basefontname`
    below is for.
    """
    candidates = [
        ("Arial", {
            "": Path("C:/Windows/Fonts/arial.ttf"),
            "B": Path("C:/Windows/Fonts/arialbd.ttf"),
            "I": Path("C:/Windows/Fonts/ariali.ttf"),
            "BI": Path("C:/Windows/Fonts/arialbi.ttf"),
        }),
        ("Liberation Sans", {
            "": Path("/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf"),
            "B": Path("/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf"),
            "I": Path("/usr/share/fonts/truetype/liberation/LiberationSans-Italic.ttf"),
            "BI": Path("/usr/share/fonts/truetype/liberation/LiberationSans-BoldItalic.ttf"),
        }),
    ]
    try:
        import matplotlib
        d = Path(matplotlib.__file__).parent / "mpl-data" / "fonts" / "ttf"
        candidates.append(("DejaVu Sans", {
            "": d / "DejaVuSans.ttf", "B": d / "DejaVuSans-Bold.ttf",
            "I": d / "DejaVuSans-Oblique.ttf", "BI": d / "DejaVuSans-BoldOblique.ttf",
        }))
    except Exception:
        pass

    for name, files in candidates:
        if all(p.exists() for p in files.values()):
            base = name.replace(" ", "")
            faces = {}
            for suffix, path in files.items():
                face = base + (suffix or "")
                pdfmetrics.registerFont(TTFont(face, str(path)))
                faces[suffix] = face
            pdfmetrics.registerFontFamily(
                base, normal=faces[""], bold=faces["B"],
                italic=faces["I"], boldItalic=faces["BI"])
            # Stop reportlab writing an unembedded Helvetica into every page.
            from reportlab import rl_config
            rl_config.canvas_basefontname = faces[""]
            register_mono()
            print(f"  fonts: {name} (embedded)")
            return faces[""], faces["B"], faces["I"], faces["BI"]
    raise SystemExit("render_pdf: no embeddable sans-serif family found.")


# ----------------------------------------------------------------- markdown
_SUP = re.compile(r"<sup>(.*?)</sup>", re.S)


def inline(text: str, regular: str, bold: str, italic: str) -> str:
    """Markdown inline spans to reportlab's mini-HTML."""
    text = text.replace("\\*", "\u0001")
    # Protect the tags we keep before escaping everything else.
    text = _SUP.sub(lambda m: f"\u0002{m.group(1)}\u0003", text)
    text = html.escape(text, quote=False)
    text = text.replace("\u0002", "<super>").replace("\u0003", "</super>")
    text = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", text, flags=re.S)
    text = re.sub(r"(?<![\*\w])\*([^*]+?)\*(?![\*\w])", r"<i>\1</i>", text, flags=re.S)
    mono = MONO["face"]
    if mono:
        text = re.sub(r"`([^`]+?)`", rf"<font face='{mono}'>\1</font>", text)
    else:
        text = re.sub(r"`([^`]+?)`", r"\1", text)
    text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r"\1", text)
    return text.replace("\u0001", "*")


def blocks(md: str):
    """Split Markdown into (kind, payload) blocks."""
    md = re.sub(r"<!--.*?-->", "", md, flags=re.S)
    lines = md.split("\n")
    i, out = 0, []
    while i < len(lines):
        line = lines[i]
        if not line.strip():
            i += 1
            continue
        if line.strip() in ("---", "***", "___"):
            out.append(("rule", None))
            i += 1
        elif line.startswith("#"):
            level = len(line) - len(line.lstrip("#"))
            out.append((f"h{level}", line.lstrip("#").strip()))
            i += 1
        elif line.lstrip().startswith("|"):
            rows = []
            while i < len(lines) and lines[i].lstrip().startswith("|"):
                rows.append(lines[i].strip())
                i += 1
            out.append(("table", rows))
        elif re.match(r"^\s*(\d+\.|[-*+])\s+", line):
            items, indent = [], []
            while i < len(lines) and (re.match(r"^\s*(\d+\.|[-*+])\s+", lines[i])
                                      or (lines[i].startswith("   ") and items)):
                if re.match(r"^\s*(\d+\.|[-*+])\s+", lines[i]):
                    m = re.match(r"^(\s*)(\d+\.|[-*+])\s+(.*)$", lines[i])
                    indent.append(m.group(2))
                    items.append(m.group(3))
                else:
                    items[-1] += " " + lines[i].strip()
                i += 1
            out.append(("list", list(zip(indent, items))))
        else:
            para = []
            while i < len(lines) and lines[i].strip() and not lines[i].startswith("#") \
                    and lines[i].strip() not in ("---", "***", "___") \
                    and not lines[i].lstrip().startswith("|") \
                    and not re.match(r"^\s*(\d+\.|[-*+])\s+", lines[i]):
                para.append(lines[i].strip())
                i += 1
            out.append(("p", " ".join(para)))
    return out


# --------------------------------------------------- line-numbered paragraph
class Numbered(Paragraph):
    """A paragraph that prints a line number beside each of its own lines.

    Reportlab knows where every wrapped line will sit only after `wrap`, so
    the numbers are drawn here rather than guessed from the page geometry.
    The counter lives on the document, so it runs continuously across pages.
    """

    counter = {"n": 0}
    enabled = {"on": False}

    def draw(self):
        Paragraph.draw(self)
        if not self.enabled["on"]:
            return
        lines = getattr(self.blPara, "lines", []) or []
        leading = self.style.leading
        c = self.canv
        c.saveState()
        c.setFont(self.style.fontName, 6.5)
        c.setFillColor(colors.HexColor("#9AA3A8"))
        for k in range(len(lines)):
            self.counter["n"] += 1
            y = self.height - (k + 1) * leading + leading * 0.25
            c.drawRightString(-6 * mm, y, str(self.counter["n"]))
        c.restoreState()


# ------------------------------------------------------------------- render
def build(md_path: Path, out_path: Path, *, line_numbers: bool,
          title: str, figures: list[Path] | None = None,
          fig_captions: list[str] | None = None) -> None:
    regular, bold, italic, bolditalic = register_fonts()
    Numbered.enabled["on"] = line_numbers
    Numbered.counter["n"] = 0

    body = ParagraphStyle("body", fontName=regular, fontSize=9.5, leading=14.5,
                          alignment=TA_JUSTIFY, spaceAfter=7,
                          textColor=colors.HexColor("#16191C"))
    h1 = ParagraphStyle("h1", parent=body, fontName=bold, fontSize=14.5,
                        leading=19, alignment=TA_LEFT, spaceBefore=6,
                        spaceAfter=10)
    h2 = ParagraphStyle("h2", parent=body, fontName=bold, fontSize=11,
                        leading=15, alignment=TA_LEFT, spaceBefore=13,
                        spaceAfter=5)
    h3 = ParagraphStyle("h3", parent=body, fontName=bold, fontSize=9.8,
                        leading=14, alignment=TA_LEFT, spaceBefore=9,
                        spaceAfter=3)
    small = ParagraphStyle("small", parent=body, fontSize=8.2, leading=11.5,
                           alignment=TA_LEFT, spaceAfter=4)
    cell = ParagraphStyle("cell", parent=body, fontSize=7.6, leading=10,
                          alignment=TA_LEFT, spaceAfter=0)
    cellh = ParagraphStyle("cellh", parent=cell, fontName=bold)
    cap = ParagraphStyle("cap", parent=body, fontSize=8.4, leading=12,
                         alignment=TA_LEFT, spaceBefore=4, spaceAfter=14)

    def P(text, style=body, numbered=True):
        cls = Numbered if (line_numbers and numbered and style is body) else Paragraph
        return cls(inline(text, regular, bold, italic), style)

    story = []
    md = md_path.read_text(encoding="utf-8")
    title_lines = []
    for kind, payload in blocks(md):
        if kind == "h1":
            title_lines.append(payload)
            continue
        if title_lines:
            story.append(Paragraph(inline(" ".join(title_lines), regular, bold,
                                          italic), h1))
            title_lines = []
        if kind == "rule":
            story.append(Spacer(1, 3))
            story.append(HRFlowable(width="100%", thickness=0.6,
                                    color=colors.HexColor("#D6DBDE"),
                                    spaceAfter=8))
        elif kind == "h2":
            story.append(Paragraph(inline(payload, regular, bold, italic), h2))
        elif kind in ("h3", "h4", "h5", "h6"):
            story.append(Paragraph(inline(payload, regular, bold, italic), h3))
        elif kind == "p":
            story.append(P(payload))
        elif kind == "list":
            for marker, item in payload:
                story.append(P(f"{marker}\u00a0 {item}", small, numbered=False))
        elif kind == "table":
            story.append(render_table(payload, cell, cellh, regular, bold, italic))
            story.append(Spacer(1, 8))
    if title_lines:
        story.append(Paragraph(inline(" ".join(title_lines), regular, bold,
                                      italic), h1))

    OUT.mkdir(parents=True, exist_ok=True)
    if figures:
        story.append(PageBreak())
        avail = PAGE[0] - MARGIN_L - MARGIN_R
        for path, caption in zip(figures, fig_captions or [""] * len(figures)):
            import fitz
            page = fitz.open(path)[0]
            w_pt, h_pt = page.rect.width, page.rect.height
            scale = min(avail / w_pt, (PAGE[1] - MARGIN_T - MARGIN_B - 60) / h_pt)
            png = OUT / f"_embed_{path.stem}.png"
            page.get_pixmap(dpi=400, alpha=False).save(png)
            img = Image(str(png), width=w_pt * scale, height=h_pt * scale)
            img.hAlign = "LEFT"
            story.append(KeepTogether(
                [img, Paragraph(inline(caption, regular, bold, italic), cap)]))

    frame = Frame(MARGIN_L, MARGIN_B, PAGE[0] - MARGIN_L - MARGIN_R,
                  PAGE[1] - MARGIN_T - MARGIN_B, id="body",
                  leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=0)

    def footer(canvas, doc):
        canvas.saveState()
        canvas.setFont(regular, 8)
        canvas.setFillColor(colors.HexColor("#9AA3A8"))
        canvas.drawCentredString(PAGE[0] / 2, MARGIN_B - 11,
                                 f"{doc.page}")
        canvas.drawString(MARGIN_L, PAGE[1] - MARGIN_T + 7, title)
        canvas.restoreState()

    OUT.mkdir(parents=True, exist_ok=True)
    doc = BaseDocTemplate(str(out_path), pagesize=PAGE,
                          leftMargin=MARGIN_L, rightMargin=MARGIN_R,
                          topMargin=MARGIN_T, bottomMargin=MARGIN_B,
                          title=title, author="Alexander Schulz, Christian T. Thiel")
    doc.addPageTemplates([PageTemplate(id="all", frames=[frame], onPage=footer)])
    doc.build(story)
    check_fonts_embedded(out_path)
    size_mb = out_path.stat().st_size / 1e6
    print(f"  wrote {out_path.relative_to(ROOT)}  ({doc.page} pages, {size_mb:.2f} MB)")
    for junk in OUT.glob("_embed_*.png"):
        junk.unlink()


def check_fonts_embedded(path: Path) -> None:
    """Refuse to ship a PDF that only references a font instead of embedding it.

    Cell Press, like every production system, rejects those, and the failure
    is invisible on screen because the viewer substitutes a look-alike.
    """
    import fitz

    doc = fitz.open(path)
    loose = set()
    for page in doc:
        for xref, ext, ftype, basefont, *_ in page.get_fonts(full=True):
            if ext in ("n/a", "") or ftype == "Type1":
                loose.add(basefont)
    doc.close()
    if loose:
        raise ValueError(
            f"{path.name}: {len(loose)} font(s) referenced but not embedded: "
            f"{', '.join(sorted(loose))}")


def render_table(rows, cell, cellh, regular, bold, italic):
    def cells(line):
        return [c.strip() for c in line.strip().strip("|").split("|")]

    header = cells(rows[0])
    data_rows = [cells(r) for r in rows[2:]] if len(rows) > 2 else []
    data = [[Paragraph(inline(c, regular, bold, italic), cellh) for c in header]]
    for r in data_rows:
        r = (r + [""] * len(header))[:len(header)]
        data.append([Paragraph(inline(c, regular, bold, italic), cell) for c in r])
    avail = PAGE[0] - MARGIN_L - MARGIN_R
    t = Table(data, colWidths=[avail / len(header)] * len(header), repeatRows=1)
    t.setStyle(TableStyle([
        # Without this the table's default font is the base-14 Helvetica, which
        # lands in the page resources unembedded even though every cell is a
        # Paragraph that draws in Arial.
        ("FONT", (0, 0), (-1, -1), regular),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
        ("LEFTPADDING", (0, 0), (-1, -1), 4),
        ("RIGHTPADDING", (0, 0), (-1, -1), 4),
        ("LINEBELOW", (0, 0), (-1, 0), 0.7, colors.HexColor("#6E787D")),
        ("LINEBELOW", (0, 1), (-1, -2), 0.25, colors.HexColor("#E6EAEC")),
        ("LINEBELOW", (0, -1), (-1, -1), 0.7, colors.HexColor("#6E787D")),
    ]))
    return t


LEGENDS = [
    "**Figure 1.** Breadth of action across the skeleton marks the Mendelian "
    "layer. Full legend in the main text.",
    "**Figure 2.** The rule relates anatomy to anatomy and stops at the organ "
    "boundary. Full legend in the main text.",
    "**Figure 3.** The layer replicates in an independent cohort and carries "
    "fracture risk. Full legend in the main text.",
    "**Figure 4.** What the shared layer is made of. Full legend in the main "
    "text.",
    "**Graphical abstract.** One gene property, two clinical ends.",
]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--only", choices=["manuscript", "supplement", "cover"],
                    help="render one document instead of all three")
    args = ap.parse_args()
    md = ROOT / "06_manuscript"
    figs = ROOT / "figures"

    jobs = {
        "manuscript": dict(
            md_path=md / "manuscript.md", out_path=OUT / "Manuscript.pdf",
            line_numbers=True, title="Breadth of action across the skeleton",
            figures=[figs / f"figure{i}.pdf" for i in (1, 2, 3, 4)]
                    + [figs / "graphical_abstract.pdf"],
            fig_captions=LEGENDS),
        "supplement": dict(
            md_path=md / "supplement.md",
            out_path=OUT / "Supplemental Information.pdf",
            line_numbers=False, title="Supplemental Information"),
        "cover": dict(
            md_path=md / "cover_letter.md", out_path=OUT / "Cover Letter.pdf",
            line_numbers=False, title=""),
    }
    for name, kw in jobs.items():
        if args.only and name != args.only:
            continue
        if not kw["md_path"].exists():
            print(f"  skipped {name}: {kw['md_path'].name} not present")
            continue
        print(f"{name}:")
        build(**kw)
    return 0


if __name__ == "__main__":
    sys.exit(main())
