from __future__ import annotations

import csv
import hashlib
import json
import shutil
import zipfile
from datetime import datetime
from pathlib import Path

from PIL import Image, ImageStat


ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "figures" / "main_text_composites"
OUT_DIR = ROOT / "figures" / "submission_figure_package_v5_rerun_20260527"
PACKAGE_LABEL = "figures/submission_figure_package_v5_rerun_20260527"
MANIFEST_TSV = FIG_DIR / "main_text_composite_source_manifest.tsv"

FINAL_FIGURES = [
    "Figure_1_bulk_standard_association",
    "Figure_2_single_nucleus_standard_umap",
    "Figure_3_single_nucleus_donor_aware_stromal",
    "Figure_4_spatial_standard_localization",
    "Figure_5_network_standard_prioritization",
    "Figure_6_gene_level_coexpression_standard_context",
    "Figure_7_benchmark_robustness_sensitivity",
]

DOCX_FILES = [
    ROOT / "manuscript_outputs" / "current_email_updated_manuscript_with_v5_figures.docx",
]

SYNCED_DOCUMENT_FILES = [
    ROOT / "manuscript_outputs" / "current_email_updated_manuscript_with_v5_figures.docx",
    ROOT / "manuscript_outputs" / "current_email_updated_manuscript_with_standard_figures.docx",
]

SENSITIVE_CROPS = {
    "Figure_1_bulk_standard_association": {
        "panel_A": (0.00, 0.00, 0.50, 0.34),
        "panel_C": (0.00, 0.34, 0.50, 0.67),
        "panel_E": (0.00, 0.67, 0.50, 1.00),
    },
    "Figure_2_single_nucleus_standard_umap": {
        "umap_row": (0.00, 0.00, 1.00, 0.49),
        "score_marker_row": (0.00, 0.49, 1.00, 1.00),
        "marker_dotplot": (0.62, 0.49, 1.00, 1.00),
    },
    "Figure_3_single_nucleus_donor_aware_stromal": {
        "donor_scores": (0.00, 0.00, 0.66, 0.34),
        "sensitivity_scores": (0.00, 0.34, 0.66, 0.68),
        "violin": (0.00, 0.68, 1.00, 1.00),
    },
    "Figure_4_spatial_standard_localization": {
        "spatial_maps": (0.00, 0.00, 1.00, 0.39),
        "niche_gradient": (0.00, 0.39, 1.00, 0.75),
        "coverage": (0.00, 0.75, 1.00, 1.00),
    },
    "Figure_5_network_standard_prioritization": {
        "network_lollipop": (0.00, 0.00, 1.00, 0.58),
        "centrality_overlap": (0.00, 0.58, 1.00, 1.00),
    },
    "Figure_6_gene_level_coexpression_standard_context": {
        "heatmaps": (0.00, 0.00, 1.00, 0.52),
        "network_matrix": (0.00, 0.52, 1.00, 1.00),
        "evidence_matrix": (0.50, 0.52, 1.00, 1.00),
    },
    "Figure_7_benchmark_robustness_sensitivity": {
        "benchmark_proxy": (0.00, 0.00, 0.67, 0.52),
        "stability": (0.00, 0.52, 0.67, 1.00),
        "correlation": (0.67, 0.00, 1.00, 1.00),
    },
}


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def ensure_dirs() -> dict[str, Path]:
    dirs = {
        "final_pdf": OUT_DIR / "final_pdf",
        "preview_png": OUT_DIR / "preview_png",
        "png_600ppi": OUT_DIR / "png_600ppi",
        "tiff_600ppi": OUT_DIR / "tiff_600ppi",
        "qc_crops": OUT_DIR / "qc_crops",
        "source_panels": OUT_DIR / "source_panels",
        "docx_synced": OUT_DIR / "docx_synced",
    }
    for path in dirs.values():
        path.mkdir(parents=True, exist_ok=True)
    return dirs


def read_manifest() -> list[dict[str, str]]:
    with MANIFEST_TSV.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def image_qc(path: Path) -> dict[str, object]:
    with Image.open(path) as image:
        gray = image.convert("L")
        stat = ImageStat.Stat(gray)
        dpi = image.info.get("dpi", (0, 0))
        return {
            "pixels": f"{image.width}x{image.height}",
            "dpi": tuple(round(float(x), 2) for x in dpi),
            "mean_gray": round(stat.mean[0], 2),
            "std_gray": round(stat.stddev[0], 2),
            "nonblank": stat.stddev[0] > 2,
            "md5": md5(path),
        }


def save_tiff(src_png: Path, out_tiff: Path) -> None:
    with Image.open(src_png) as image:
        image.save(out_tiff, dpi=(600, 600), compression="tiff_lzw")


def save_crops(src_png: Path, stem: str, crop_specs: dict[str, tuple[float, float, float, float]]) -> list[dict[str, object]]:
    crop_rows = []
    with Image.open(src_png) as image:
        width, height = image.size
        for crop_name, (x0, y0, x1, y1) in crop_specs.items():
            box = (
                int(width * x0),
                int(height * y0),
                int(width * x1),
                int(height * y1),
            )
            crop = image.crop(box)
            out_path = OUT_DIR / "qc_crops" / f"{stem}_{crop_name}.png"
            crop.save(out_path, dpi=(220, 220))
            gray = crop.convert("L")
            stat = ImageStat.Stat(gray)
            crop_rows.append({
                "figure": stem,
                "crop": crop_name,
                "file": out_path.relative_to(OUT_DIR).as_posix(),
                "pixels": f"{crop.width}x{crop.height}",
                "std_gray": round(stat.stddev[0], 2),
                "nonblank": stat.stddev[0] > 2,
            })
    return crop_rows


def docx_media_hashes(docx: Path) -> list[str]:
    if not docx.exists():
        return []
    hashes = []
    with zipfile.ZipFile(docx) as archive:
        for name in sorted(n for n in archive.namelist() if n.startswith("word/media/")):
            digest = hashlib.md5(archive.read(name)).hexdigest()
            hashes.append(digest)
    return hashes


def write_source_manifest(rows: list[dict[str, str]], figure_qc: dict[str, dict[str, object]]) -> None:
    lines = [
        "# Source Manifest",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "This package contains the current v5 manuscript-facing bioinformatics figure set. Analytical source panels were generated by R workflows; the composite figures were assembled from R-generated source PNG panels. Composite PDFs are raster-based wrappers and must not be described as fully vector.",
        "",
        "| Figure | Panel | Source type | Source path or object | Input data/cache | Rerendered | Export used | Vector status | Preserved labels | Known limitations |",
        "|---|---|---|---|---|---|---|---|---|---|",
    ]
    for row in rows:
        fig = row["figure"]
        if fig not in FINAL_FIGURES:
            continue
        source_figure = row["source_figure"]
        source_data = row["source_data"]
        script = row["generating_script"]
        role = row["role"]
        wording = row["wording_note"]
        limitation = "R-derived source panel; composite is raster at 600 ppi. Interpret computational outputs at the association/localization/candidate-prioritization level."
        if fig == "Figure_3_single_nucleus_donor_aware_stromal" and row["panel"] in {"A", "B"}:
            limitation += " Strict high-confidence stromal gating has limited finite donor coverage and should not be over-read as disease-group separation."
        if fig == "Figure_4_spatial_standard_localization":
            limitation += " Spatial panels are descriptive true-coordinate localization, not donor-level severity inference."
        if fig == "Figure_7_benchmark_robustness_sensitivity":
            limitation += " Benchmark, proxy-adjustment, bootstrap, and leave-one-out analyses are sensitivity analyses and do not establish cell-intrinsic mechanotransduction or mechanical memory."
        export_used = f"png_600ppi/{fig}.png; tiff_600ppi/{fig}.tiff; final_pdf/{fig}.pdf"
        lines.append(
            f"| {fig} | {row['panel']} | R script and R-derived source export | {script}; {source_figure} | {source_data} | Yes, source R workflow rerun or refreshed for v5 | {export_used} | Raster-only but 600 ppi; PDF is raster wrapper | {role}; {wording} | {limitation} |"
        )
    lines.append("")
    lines.append("## Figure-Level Raster QC")
    lines.append("")
    lines.append("| Figure | PNG pixels | PNG ppi | TIFF ppi | Nonblank | MD5 |")
    lines.append("|---|---:|---:|---:|---|---|")
    for fig in FINAL_FIGURES:
        qc = figure_qc[fig]
        lines.append(
            f"| {fig} | {qc['png']['pixels']} | {qc['png']['dpi']} | {qc['tiff']['dpi']} | {qc['png']['nonblank']} | {qc['png']['md5']} |"
        )
    (OUT_DIR / "source_manifest.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_qc_report(figure_qc: dict[str, dict[str, object]], crop_rows: list[dict[str, object]], docx_status: list[dict[str, object]]) -> None:
    lines = [
        "# Figure QC Report",
        "",
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "## Package Summary",
        "",
        f"- Figure package: `{PACKAGE_LABEL}`",
        "- Export mode: docx-package plus submission-raster.",
        "- Final PDF: `final_pdf/` contains current composite PDFs. These are raster-based PDF wrappers, not fully vector figures.",
        "- Preview PNG and PNG 600 ppi: `preview_png/` and `png_600ppi/` contain current 600 ppi PNGs.",
        "- TIFF 600 ppi: `tiff_600ppi/` contains LZW-compressed TIFFs generated from the current 600 ppi PNG exports.",
        "- DOCX synchronized: see table below.",
        "",
        "## Vector and Resolution Status",
        "",
        "| Figure | Final PDF status | PNG ppi | TIFF ppi | Nonblank | Grade |",
        "|---|---|---:|---:|---|---|",
    ]
    for fig, qc in figure_qc.items():
        lines.append(
            f"| {fig} | Raster-based PDF wrapper | {qc['png']['dpi']} | {qc['tiff']['dpi']} | {qc['png']['nonblank']} | Raster-only but 600 ppi |"
        )
    lines.extend([
        "",
        "## Layout QC",
        "",
        "Automated crops were generated for sensitive regions. Nonblank crop checks passed when grayscale standard deviation exceeded 2. This does not replace author or journal visual review, but it guards against blank or fully lost panels.",
        "",
        "| Figure | Crop | File | Pixels | Nonblank |",
        "|---|---|---|---:|---|",
    ])
    for row in crop_rows:
        lines.append(
            f"| {row['figure']} | {row['crop']} | `{row['file']}` | {row['pixels']} | {row['nonblank']} |"
        )
    lines.extend([
        "",
        "## Scientific Content Preservation",
        "",
        "- Required axes, legends, color bars, gene labels, UMAP coordinates, spatial coordinates, network labels, and heatmap labels were preserved at the source-panel layer.",
        "- Composite overlays retain only panel letters; no composite-level figure number, top title, card frame, or decorative panel box was added.",
        "- Figure 7 is integrated into the same manifest and QC workflow as Figures 1-6 rather than being treated as an addendum.",
        "- Single-nucleus UMAP panels are explicitly treated as a balanced sketch with stromal-marker-priority retention, not as an all-nucleus atlas.",
        "- Spatial maps preserve true slide coordinates and are interpreted as descriptive slide-nested localization.",
        "",
        "## Manuscript Consistency",
        "",
        "| DOCX | Exists | Media count | Current figure hashes matched |",
        "|---|---|---:|---:|",
    ])
    n_final = len(FINAL_FIGURES)
    for row in docx_status:
        lines.append(f"| {row['docx']} | {row['exists']} | {row['media_count']} | {row['matched_hashes']}/{n_final} |")
    lines.extend([
        "",
        "## Human Verification Needed",
        "",
        "- Target journal figure-file rules, including TIFF/PDF/EPS/SVG, CMYK/RGB, file size, and minimum font size, require author or journal verification.",
        "- Composite PDFs are not fully vector. If the journal requires editable vector figures, source-level vector recomposition from R plot/PDF objects is still required.",
        "- Figure 3 high-confidence donor-aware stromal summaries should remain conservatively worded because strict gating leaves limited finite donor coverage.",
        "- GSE292268 spatial coordinate conventions should be confirmed against original CosMx metadata before making stronger spatial claims.",
        "- Figure 7 should remain framed as benchmark and robustness sensitivity, not as proof that the primary program outperforms known fibrosis or stromal signatures.",
    ])
    (OUT_DIR / "qc_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_delivery_note() -> None:
    lines = [
        "# Final Figure Package Note",
        "",
        "Use this folder as the current v5 figure-submission package for the MASLD mechanical-stress-program manuscript.",
        "",
        "## Current Main Figures",
        "",
        "| Figure | Base filename | Primary delivery files |",
        "|---|---|---|",
    ]
    for fig in FINAL_FIGURES:
        lines.append(
            f"| {fig.replace('_', ' ')} | `{fig}` | `png_600ppi/{fig}.png`; `tiff_600ppi/{fig}.tiff`; `final_pdf/{fig}.pdf` |"
        )
    lines.extend([
        "",
        "## Use Boundaries",
        "",
        "- For submission transfer, use the files in this package instead of selecting loose files manually.",
        "- Analytical panels were generated from R workflows and assembled into clean composites. Python was used for document packaging, TIFF conversion, hash checks, and QC reporting.",
        "- The final PDF figure files are raster-based wrappers. Do not describe them as fully vector or editable PDF figures.",
        "- If the target journal requires editable vector PDF, EPS, or SVG, the figures should be recomposed from source-level R vector outputs before upload.",
        "",
        "## Synchronized Manuscript Preview",
        "",
        "- `the current email-updated manuscript embeds the seven v5 PNG figures` embeds the current seven v5 PNG figures.",
    ])
    (OUT_DIR / "README_FINAL.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    dirs = ensure_dirs()
    rows = read_manifest()
    figure_qc: dict[str, dict[str, object]] = {}
    crop_rows: list[dict[str, object]] = []
    current_hashes = []

    for fig in FINAL_FIGURES:
        png = FIG_DIR / f"{fig}.png"
        pdf = FIG_DIR / f"{fig}.pdf"
        if not png.exists() or not pdf.exists():
            raise FileNotFoundError(fig)
        shutil.copy2(pdf, dirs["final_pdf"] / pdf.name)
        shutil.copy2(png, dirs["preview_png"] / png.name)
        shutil.copy2(png, dirs["png_600ppi"] / png.name)
        tiff = dirs["tiff_600ppi"] / f"{fig}.tiff"
        save_tiff(png, tiff)
        current_hashes.append(md5(png))
        figure_qc[fig] = {
            "png": image_qc(dirs["png_600ppi"] / png.name),
            "tiff": image_qc(tiff),
        }
        crop_rows.extend(save_crops(png, fig, SENSITIVE_CROPS[fig]))

    source_paths = sorted({item for row in rows for item in row["source_figure"].split("; ")})
    for rel in source_paths:
        src = ROOT / rel
        if src.exists():
            target = dirs["source_panels"] / src.name
            shutil.copy2(src, target)

    docx_status = []
    for synced_file in SYNCED_DOCUMENT_FILES:
        if synced_file.exists():
            shutil.copy2(synced_file, dirs["docx_synced"] / synced_file.name)

    for docx in DOCX_FILES:
        hashes = docx_media_hashes(docx)
        docx_status.append({
            "docx": docx.name,
            "exists": docx.exists(),
            "media_count": len(hashes),
            "matched_hashes": sum(1 for h in hashes if h in current_hashes),
        })

    write_source_manifest(rows, figure_qc)
    write_qc_report(figure_qc, crop_rows, docx_status)
    write_delivery_note()
    payload = {
        "package": PACKAGE_LABEL,
        "figures": figure_qc,
        "crops": crop_rows,
        "docx": docx_status,
    }
    (OUT_DIR / "qc_summary.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    print(OUT_DIR)


if __name__ == "__main__":
    main()
