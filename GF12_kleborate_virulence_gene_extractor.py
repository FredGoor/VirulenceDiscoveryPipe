#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Kleborate-guided Klebsiella virulence gene extractor

Purpose
-------
Run Kleborate on a Klebsiella pneumoniae species complex genome assembly and
extract strain-specific locus tags for virulence-associated gene families from a
matching GFF/GFF3 or GenBank annotation file.

Kleborate is primarily a locus/module typer. This script therefore uses
Kleborate to detect Klebsiella-specific virulence modules, then searches the
provided annotation file to export the corresponding genes and locus tags.

Extracted categories
--------------------
1. Capsule / K locus / surface polysaccharide
2. LPS / O-antigen
3. Yersiniabactin
4. Aerobactin
5. Salmochelin
6. Enterobactin / ferrienterobactin
7. Colibactin / genotoxin
8. Hypermucoidy / rmp regulators
9. Fimbriae / pili / adhesins
10. T6SS
11. Flagella / motility / chemotaxis
12. Toxins / toxin-antitoxin / competition toxins
13. Other surface / immune evasion candidates

Outputs
-------
The Excel workbook contains:

    - Kleborate_raw
    - Kleborate_modules
    - Extracted_genes
    - Summary_counts
    - Locus_tag_lists
    - Notes
    - one sheet per extracted category

Requirements
------------
Install the command-line and Python dependencies in the active environment:

    conda install -c conda-forge pandas openpyxl biopython
    # Install Kleborate following the instructions for your platform/environment.

Recommended use
---------------
Run interactively with file pickers:

    conda activate virulence_screen
    python GF12_kleborate_virulence_gene_extractor_publishable.py

Or provide paths directly:

    python GF12_kleborate_virulence_gene_extractor_publishable.py \
        --assembly genome.fna \
        --annotation annotation.gff \
        --output-dir results
"""


# =============================================================================
# USER SETTINGS
# =============================================================================

RUN_KLEBORATE = True

# Kleborate command. Usually simply "kleborate" if installed in the active conda env.
KLEBORATE_EXE = "kleborate"

# Kleborate preset for Klebsiella pneumoniae species complex.
KLEBORATE_PRESET = "kpsc"

# If True, only extract categories supported by a positive Kleborate module call
# for ybt/clb/iuc/iro/rmp. Categories such as fimbriae, T6SS, flagella, toxins,
# capsule and O-antigen are still extracted because Kleborate does not output
# per-gene lists for all of these categories.
REQUIRE_KLEBORATE_MODULE_FOR_ACQUIRED_VIRULENCE_LOCI = True

# Keep borderline/general hits. Useful for discovery, but can be set to False
# for a stricter list.
INCLUDE_BORDERLINE = True

# Minimum overlap is not used here because extraction is annotation-pattern-based.
# Locus tags are taken from the current locus_tag field in the annotation file.
PREFER_CURRENT_LOCUS_TAG = True

# Output workbook name
OUTPUT_EXCEL_NAME = "Kleborate_extracted_Klebsiella_virulence_genes.xlsx"

# Output command log
KLEBORATE_LOG_NAME = "kleborate_run.log"

# =============================================================================
# IMPORTS
# =============================================================================

import os
import re
import sys
import json
import shlex
import shutil
import argparse
import subprocess
from pathlib import Path
from urllib.parse import unquote

import pandas as pd

try:
    import tkinter as tk
    from tkinter import filedialog, messagebox
except Exception:  # Tkinter may be unavailable on headless systems.
    tk = None
    filedialog = None
    messagebox = None


# =============================================================================
# FILE PICKERS
# =============================================================================

def _tk_available():
    return tk is not None and filedialog is not None


def choose_file(title, patterns):
    """Open a file picker and return the selected path."""
    if not _tk_available():
        raise RuntimeError(
            "Tkinter file picker is unavailable. Provide paths on the command line "
            "with --assembly, --annotation and --output-dir."
        )

    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename(
        title=title,
        filetypes=patterns,
    )
    root.destroy()

    if not file_path:
        raise SystemExit(f"No file selected for: {title}")
    return Path(file_path)


def choose_directory(title):
    """Open a directory picker and return the selected path."""
    if not _tk_available():
        raise RuntimeError(
            "Tkinter directory picker is unavailable. Provide --output-dir on the command line."
        )

    root = tk.Tk()
    root.withdraw()
    dir_path = filedialog.askdirectory(title=title)
    root.destroy()

    if not dir_path:
        raise SystemExit(f"No directory selected for: {title}")
    return Path(dir_path)


def parse_command_line_args():
    parser = argparse.ArgumentParser(
        description="Kleborate-guided extraction of Klebsiella virulence-associated locus tags."
    )
    parser.add_argument(
        "--assembly",
        default="",
        help="Input Klebsiella genome assembly FASTA file.",
    )
    parser.add_argument(
        "--annotation",
        default="",
        help="Matching GFF/GFF3/GBK/GBFF annotation file.",
    )
    parser.add_argument(
        "--output-dir",
        default="",
        help="Output directory.",
    )
    parser.add_argument(
        "--no-gui",
        action="store_true",
        help="Disable file pickers and require all paths as command-line arguments.",
    )
    return parser.parse_args()


def resolve_runtime_paths():
    """Resolve assembly, annotation and output paths from CLI arguments or file pickers."""
    args = parse_command_line_args()

    assembly_path = str(args.assembly).strip()
    annotation_path = str(args.annotation).strip()
    outdir = str(args.output_dir).strip()

    if not assembly_path:
        if args.no_gui:
            raise SystemExit("Missing --assembly.")
        assembly_path = str(choose_file(
            "Select Klebsiella genome assembly FASTA",
            [("FASTA files", "*.fna *.fa *.fasta *.fna.gz *.fa.gz *.fasta.gz"), ("All files", "*.*")]
        ))

    if not annotation_path:
        if args.no_gui:
            raise SystemExit("Missing --annotation.")
        annotation_path = str(choose_file(
            "Select matching annotation file",
            [("Annotation files", "*.gff *.gff3 *.gbk *.gbff *.gb"), ("All files", "*.*")]
        ))

    if not outdir:
        if args.no_gui:
            raise SystemExit("Missing --output-dir.")
        outdir = str(choose_directory("Select output folder"))

    return Path(assembly_path), Path(annotation_path), Path(outdir)


# =============================================================================
# PATH / TEXT HELPERS
# =============================================================================

def clean_string(x):
    if x is None or pd.isna(x):
        return ""
    return str(x).strip()


def normalize_text(x):
    s = clean_string(x).lower()
    s = s.replace("_", " ")
    s = s.replace("-", " ")
    s = re.sub(r"[^a-z0-9]+", " ", s)
    s = re.sub(r"\s+", " ", s).strip()
    return s


def safe_sheet_name(s):
    s = re.sub(r"[^A-Za-z0-9]+", "_", str(s)).strip("_")
    return s[:31] if s else "Sheet"


def safe_sample_name(path):
    name = Path(path).stem
    name = re.sub(r"[^A-Za-z0-9_.-]+", "_", name)
    return name.strip("_")


def is_blank_value(x):
    s = clean_string(x)
    return s == "" or s.lower() in {"nan", "none", "na", "n/a", "-", "not detected"}


def command_exists(cmd):
    return shutil.which(cmd) is not None


# =============================================================================
# ANNOTATION PARSING
# =============================================================================

def parse_gff_attributes(attr_text):
    """
    Parse GFF3/GFF-like attributes.
    Handles key=value and key "value" styles where possible.
    """
    attrs = {}
    text = clean_string(attr_text)

    for part in text.split(";"):
        part = part.strip()
        if not part:
            continue

        if "=" in part:
            key, value = part.split("=", 1)
        elif " " in part:
            key, value = part.split(" ", 1)
            value = value.strip().strip('"')
        else:
            continue

        key = key.strip()
        value = unquote(value.strip().strip('"'))
        attrs[key] = value

    return attrs


def parse_gff_annotation(gff_path):
    """
    Parse CDS/gene features from GFF/GFF3 and return a dataframe.

    Output columns:
        contig, start, end, strand, feature_type, locus_tag, old_locus_tag,
        gene, product, note, raw_attributes
    """
    rows = []

    with open(gff_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            contig, source, feature_type, start, end, score, strand, phase, attrs_raw = parts

            if feature_type.lower() not in {"cds", "gene"}:
                continue

            attrs = parse_gff_attributes(attrs_raw)

            locus_tag = (
                attrs.get("locus_tag")
                or attrs.get("locus-tag")
                or attrs.get("Locus_tag")
                or attrs.get("ID")
                or attrs.get("Name")
                or ""
            )

            old_locus_tag = (
                attrs.get("old_locus_tag")
                or attrs.get("old-locus-tag")
                or attrs.get("Old_locus_tag")
                or ""
            )

            gene = attrs.get("gene") or attrs.get("Name") or ""
            product = attrs.get("product") or attrs.get("Product") or ""
            note = attrs.get("note") or attrs.get("Note") or ""

            rows.append({
                "contig": contig,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "feature_type": feature_type,
                "locus_tag": locus_tag,
                "old_locus_tag": old_locus_tag,
                "gene": gene,
                "product": product,
                "note": note,
                "raw_attributes": attrs_raw,
            })

    df = pd.DataFrame(rows)

    if df.empty:
        raise ValueError(f"No CDS/gene features parsed from GFF: {gff_path}")

    return df


def parse_genbank_annotation(gbk_path):
    """
    Parse CDS/gene features from GenBank/GBFF using Biopython.

    If Biopython is not installed:
        conda install -c conda-forge biopython
    """
    try:
        from Bio import SeqIO
    except ImportError as exc:
        raise ImportError(
            "Biopython is required to parse .gbk/.gbff files. "
            "Install it with: conda install -c conda-forge biopython"
        ) from exc

    rows = []

    for record in SeqIO.parse(str(gbk_path), "genbank"):
        contig = record.id

        for feature in record.features:
            if feature.type.lower() not in {"cds", "gene"}:
                continue

            qualifiers = feature.qualifiers
            start = int(feature.location.start) + 1
            end = int(feature.location.end)
            strand = "+" if feature.location.strand == 1 else "-" if feature.location.strand == -1 else "."

            def q1(key):
                vals = qualifiers.get(key, [])
                return vals[0] if vals else ""

            locus_tag = q1("locus_tag") or q1("protein_id") or q1("gene")
            old_locus_tag = q1("old_locus_tag")
            gene = q1("gene")
            product = q1("product")
            note = "; ".join(qualifiers.get("note", []))

            rows.append({
                "contig": contig,
                "start": start,
                "end": end,
                "strand": strand,
                "feature_type": feature.type,
                "locus_tag": locus_tag,
                "old_locus_tag": old_locus_tag,
                "gene": gene,
                "product": product,
                "note": note,
                "raw_attributes": json.dumps(qualifiers, ensure_ascii=False),
            })

    df = pd.DataFrame(rows)

    if df.empty:
        raise ValueError(f"No CDS/gene features parsed from GenBank file: {gbk_path}")

    return df


def parse_annotation(annotation_path):
    suffix = annotation_path.suffix.lower()

    if suffix in {".gff", ".gff3"}:
        return parse_gff_annotation(annotation_path)

    if suffix in {".gbk", ".gbff", ".gb"}:
        return parse_genbank_annotation(annotation_path)

    raise ValueError(
        f"Unsupported annotation format: {annotation_path}\n"
        "Use .gff, .gff3, .gbk or .gbff."
    )


# =============================================================================
# KLEBORATE
# =============================================================================

def run_kleborate(assembly_path, outdir, sample_id):
    """
    Run Kleborate and return output path plus log path.
    """
    if not command_exists(KLEBORATE_EXE):
        raise FileNotFoundError(
            f"Kleborate executable not found: {KLEBORATE_EXE}\n"
            "Activate your conda environment first, e.g. conda activate virulence_screen"
        )

    kleb_out = outdir / f"{sample_id}_kleborate_results"
    log_path = outdir / KLEBORATE_LOG_NAME

    cmd = [
        KLEBORATE_EXE,
        "-a", str(assembly_path),
        "-o", str(kleb_out),
        "-p", KLEBORATE_PRESET,
        "--trim_headers",
    ]

    print("[INFO] Running Kleborate:")
    print(" ".join(shlex.quote(c) for c in cmd))

    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    with open(log_path, "w", encoding="utf-8", errors="replace") as f:
        f.write("[COMMAND]\n")
        f.write(" ".join(shlex.quote(c) for c in cmd))
        f.write("\n\n[RETURNCODE]\n")
        f.write(str(result.returncode))
        f.write("\n\n[STDOUT]\n")
        f.write(result.stdout or "")
        f.write("\n\n[STDERR]\n")
        f.write(result.stderr or "")

    if result.returncode != 0:
        raise RuntimeError(
            f"Kleborate failed. See log:\n{log_path}\n\n"
            f"STDERR:\n{result.stderr[:2000]}"
        )

    return kleb_out, log_path


def find_kleborate_result_tables(kleb_out):
    """
    Kleborate versions may write either a file or a directory.
    Return candidate tabular result files.
    """
    kleb_out = Path(kleb_out)
    candidates = []

    if kleb_out.is_file():
        candidates.append(kleb_out)

    if kleb_out.is_dir():
        for ext in ("*.tsv", "*.txt", "*.csv"):
            candidates.extend(kleb_out.rglob(ext))

    # Also check paths with common suffixes
    for suffix in [".tsv", ".txt", ".csv"]:
        p = Path(str(kleb_out) + suffix)
        if p.exists() and p.is_file():
            candidates.append(p)

    # Keep non-empty files that look tabular
    final = []
    for p in candidates:
        if p.exists() and p.stat().st_size > 0:
            try:
                first = p.read_text(encoding="utf-8", errors="replace").splitlines()[0]
                if "\t" in first or "," in first:
                    final.append(p)
            except Exception:
                pass

    # Deduplicate
    seen = set()
    unique = []
    for p in final:
        rp = str(p.resolve())
        if rp not in seen:
            seen.add(rp)
            unique.append(p)

    return unique


def read_kleborate_table(kleb_out):
    """
    Read the most plausible Kleborate summary table.
    """
    tables = find_kleborate_result_tables(kleb_out)

    if not tables:
        raise FileNotFoundError(f"No Kleborate tabular output found at: {kleb_out}")

    parsed = []
    for p in tables:
        try:
            sep = "\t"
            df = pd.read_csv(p, sep=sep, dtype=str).fillna("")
            if df.shape[1] < 2:
                df = pd.read_csv(p, sep=",", dtype=str).fillna("")
            if df.shape[1] >= 2:
                parsed.append((p, df))
        except Exception:
            continue

    if not parsed:
        raise ValueError(f"Could not parse any Kleborate result table from: {kleb_out}")

    # Prefer table with most columns, usually the main summary.
    parsed.sort(key=lambda x: x[1].shape[1], reverse=True)
    chosen_path, chosen_df = parsed[0]

    print(f"[INFO] Parsed Kleborate table: {chosen_path}")
    return chosen_df, chosen_path


def kleborate_module_value(row, keyword_patterns):
    """
    Search Kleborate output columns for relevant module calls.
    Returns concatenated non-empty values from matching columns.
    """
    values = []
    for col in row.index:
        norm_col = normalize_text(col)
        if any(re.search(p, norm_col, flags=re.IGNORECASE) for p in keyword_patterns):
            val = clean_string(row[col])
            if not is_blank_value(val):
                values.append(f"{col}={val}")
    return "; ".join(values)


def module_is_detected(value_text):
    """
    Interpret Kleborate calls.
    Negative calls are usually 0, -, empty, none, not detected.
    Positive calls include alleles, STs, lineages, K/O types, incomplete calls, etc.
    """
    s = clean_string(value_text)
    if not s:
        return False

    # Remove column names and punctuation for a broad negative check
    lower = s.lower()
    negative_tokens = {
        "0", "-", "none", "nan", "not detected", "negative", "unknown",
        "missing", "no"
    }

    pieces = re.split(r"[;,\s]+", lower)
    nontrivial = []
    for p in pieces:
        p = p.strip()
        if not p:
            continue
        if "=" in p:
            p = p.split("=", 1)[1]
        if p not in negative_tokens:
            nontrivial.append(p)

    return len(nontrivial) > 0


def detect_kleborate_modules(kleborate_df):
    """
    Detect major virulence modules from Kleborate output columns.
    Works flexibly across Kleborate versions by matching column names.
    """
    if kleborate_df.empty:
        return pd.DataFrame()

    # Usually one genome row; if multiple, use first but keep all rows in raw output.
    row = kleborate_df.iloc[0]

    module_patterns = {
        "K locus / capsule typing": [r"\bk\b", r"k locus", r"kaptive", r"capsule", r"wzi"],
        "O locus / LPS typing": [r"\bo\b", r"o locus", r"o antigen", r"lps"],
        "Yersiniabactin": [r"ybt", r"yersiniabactin", r"ybst", r"icekp"],
        "Colibactin": [r"clb", r"colibactin", r"cbst"],
        "Aerobactin": [r"iuc", r"iut", r"aerobactin"],
        "Salmochelin": [r"iro", r"salmochelin"],
        "Hypermucoidy / rmp": [r"rmp", r"hypermuc", r"mucovisc"],
        "Virulence score": [r"virulence score", r"virulence_score"],
    }

    rows = []
    for module, patterns in module_patterns.items():
        value = kleborate_module_value(row, patterns)
        detected = module_is_detected(value)
        rows.append({
            "module": module,
            "detected_or_reported": detected,
            "kleborate_values": value,
        })

    return pd.DataFrame(rows)


# =============================================================================
# EXTRACTION RULES
# =============================================================================

CATEGORY_RULES = {
    "Capsule / K locus / surface polysaccharide": {
        "always_extract": True,
        "kleborate_module": "K locus / capsule typing",
        "core": [
            r"\bwzi\b", r"\bwza\b", r"\bwzb\b", r"\bwzc\b",
            r"\bcps\b", r"capsule", r"capsular",
            r"capsular polysaccharide", r"K locus",
            r"\bgalF\b", r"\bugd\b",
            r"capsule assembly", r"surface polysaccharide",
            r"polysaccharide export",
            r"glycosyltransferase.*capsule",
            r"ctr capsule",
        ],
        "borderline": [
            r"glycosyltransferase",
            r"galactosyl transferase",
            r"UDP.*galactose",
            r"UDP.*glucose",
            r"carbohydrate lyase",
        ],
    },

    "LPS / O-antigen": {
        "always_extract": True,
        "kleborate_module": "O locus / LPS typing",
        "core": [
            r"lipopolysaccharide", r"\blps\b",
            r"O-antigen", r"O antigen",
            r"\bwzx\b", r"\bwzy\b", r"\bwzm\b", r"\bwzt\b",
            r"\bwaa[A-Za-z0-9]*\b", r"\brfa[A-Za-z0-9]*\b", r"\brfb[A-Za-z0-9]*\b",
            r"\brml[A-Za-z0-9]*\b",
            r"heptosyltransferase", r"phosphoheptose",
            r"O-antigen ABC transport",
            r"O antigen ABC transport",
        ],
        "borderline": [
            r"glycosyltransferase",
            r"sugar transferase",
            r"polysaccharide",
        ],
    },

    "Yersiniabactin": {
        "always_extract": False,
        "kleborate_module": "Yersiniabactin",
        "core": [
            r"\bybt[A-Za-z0-9]*\b", r"yersiniabactin",
            r"\birp[0-9A-Za-z]*\b", r"\bfyuA\b",
            r"HMWP1", r"HMWP2",
            r"pesticin/yersiniabactin",
            r"salicylate synthase",
        ],
        "borderline": [
            r"TonB-dependent receptor.*yersiniabactin",
        ],
    },

    "Aerobactin": {
        "always_extract": False,
        "kleborate_module": "Aerobactin",
        "core": [
            r"\biuc[A-Za-z0-9]*\b", r"\biutA\b", r"aerobactin",
        ],
        "borderline": [],
    },

    "Salmochelin": {
        "always_extract": False,
        "kleborate_module": "Salmochelin",
        "core": [
            r"\biro[A-Za-z0-9]*\b", r"salmochelin",
        ],
        "borderline": [],
    },

    "Enterobactin / ferrienterobactin": {
        "always_extract": True,
        "kleborate_module": "",
        "core": [
            r"\bent[A-Za-z0-9]*\b", r"\bfep[A-Za-z0-9]*\b", r"\bfes\b",
            r"enterobactin", r"enterochelin", r"ferrienterobactin",
            r"ferric enterobactin", r"iron-enterobactin",
            r"isochorismate synthase",
            r"2,3-dihydro-2,3-dihydroxybenzoate",
        ],
        "borderline": [],
    },

    "Colibactin / genotoxin": {
        "always_extract": False,
        "kleborate_module": "Colibactin",
        "core": [
            r"\bclb[A-Za-z0-9]*\b", r"colibactin",
            r"genotoxin",
        ],
        "borderline": [],
    },

    "Hypermucoidy / rmp regulators": {
        "always_extract": False,
        "kleborate_module": "Hypermucoidy / rmp",
        "core": [
            r"\brmp[A-Za-z0-9]*\b", r"\brcsB\b",
            r"hypermuc", r"mucovisc",
            r"regulator.*capsule",
            r"positive regulator.*capsule",
        ],
        "borderline": [
            r"\brcs[A-Za-z0-9]*\b",
        ],
    },

    "Fimbriae / pili / adhesins": {
        "always_extract": True,
        "kleborate_module": "",
        "core": [
            r"\bmrk[A-Za-z0-9]*\b", r"\bfim[A-Za-z0-9]*\b", r"\becp[A-Za-z0-9]*\b",
            r"fimbria", r"fimbrial", r"fimbriae",
            r"pilus", r"\bpili\b", r"pilin",
            r"usher", r"fimbrial.*chaperone", r"periplasmic chaperone.*fim",
            r"adhesin", r"hemagglutinin", r"haemagglutinin",
            r"autotransporter",
            r"curli", r"\bcsg[A-Za-z0-9]*\b",
        ],
        "borderline": [
            r"putative receptor",
            r"outer membrane protein A", r"\bompA\b",
        ],
    },

    "T6SS": {
        "always_extract": True,
        "kleborate_module": "",
        "core": [
            r"type VI secretion", r"type 6 secretion", r"\bT6SS\b",
            r"\btss[A-Za-z0-9]*\b", r"\bhcp[A-Za-z0-9]*\b",
            r"\bvgr[A-Za-z0-9]*\b", r"\bvgrG\b",
            r"\bclpV\b", r"\bdotU\b", r"\bicmF\b",
            r"\bevpB\b", r"\brhs\b",
            r"contractile sheath", r"baseplate",
        ],
        "borderline": [
            r"OmpA/MotB domain",
        ],
    },

    "Flagella / motility / chemotaxis": {
        "always_extract": True,
        "kleborate_module": "",
        "core": [
            r"flagell", r"\bfli[A-Za-z0-9]*\b", r"\bflg[A-Za-z0-9]*\b",
            r"\bflh[A-Za-z0-9]*\b", r"\bmot[A-Za-z0-9]*\b",
            r"\bche[A-Za-z0-9]*\b", r"chemotaxis", r"motility",
        ],
        "borderline": [
            r"phase switching",
        ],
    },

    "Toxins / toxin-antitoxin / competition toxins": {
        "always_extract": True,
        "kleborate_module": "",
        "core": [
            r"toxin", r"antitoxin", r"toxin-antitoxin",
            r"\bRelE\b", r"\bParE\b", r"\bYafO\b", r"\bSymE\b",
            r"\bTisB\b", r"\bAbiEi\b", r"\bAbiEii\b", r"\bAbiGii\b",
            r"\bPhd\b", r"\bYefM\b",
            r"small toxic", r"toxic membrane polypeptide",
            r"host cell division inhibitor",
            r"\bCdi\b", r"Icd-like",
            r"bacteriocin", r"microcin",
        ],
        "borderline": [
            r"nuclease inhibitor",
        ],
    },

    "Other surface / immune evasion candidates": {
        "always_extract": True,
        "kleborate_module": "",
        "core": [
            r"serum resistance", r"complement resistance",
            r"\btraT\b", r"\biss\b",
            r"outer membrane protein A", r"\bompA\b",
            r"biofilm", r"\bbap\b",
        ],
        "borderline": [
            r"\bNlpI\b", r"lipoprotein NlpI",
            r"outer membrane",
        ],
    },
}


def category_allowed_by_kleborate(category, detected_modules_df):
    rules = CATEGORY_RULES[category]
    module = rules.get("kleborate_module", "")

    if rules.get("always_extract", True):
        return True, ""

    if not REQUIRE_KLEBORATE_MODULE_FOR_ACQUIRED_VIRULENCE_LOCI:
        return True, "module requirement disabled"

    if not module:
        return True, ""

    if detected_modules_df.empty:
        return False, f"Kleborate module '{module}' not available"

    hits = detected_modules_df[detected_modules_df["module"] == module]
    if hits.empty:
        return False, f"Kleborate module '{module}' not found in parsed output"

    detected = bool(hits.iloc[0]["detected_or_reported"])
    reason = hits.iloc[0]["kleborate_values"]

    return detected, reason


def match_patterns(text, patterns):
    for pat in patterns:
        if re.search(pat, text, flags=re.IGNORECASE):
            return pat
    return ""


def extract_genes_by_category(annotation_df, detected_modules_df):
    """
    Extract genes from annotation according to category rules.
    """
    rows = []

    for _, feature in annotation_df.iterrows():
        locus_tag = clean_string(feature.get("locus_tag", ""))
        old_locus_tag = clean_string(feature.get("old_locus_tag", ""))
        gene = clean_string(feature.get("gene", ""))
        product = clean_string(feature.get("product", ""))
        note = clean_string(feature.get("note", ""))
        raw = clean_string(feature.get("raw_attributes", ""))

        text = " ".join([locus_tag, old_locus_tag, gene, product, note, raw])

        if not locus_tag:
            # Avoid making up locus tags from old tags/protein IDs.
            continue

        for category, rule_set in CATEGORY_RULES.items():
            allowed, kleb_reason = category_allowed_by_kleborate(category, detected_modules_df)
            if not allowed:
                continue

            matched_rule = match_patterns(text, rule_set.get("core", []))
            confidence = "core"

            if not matched_rule and INCLUDE_BORDERLINE:
                matched_rule = match_patterns(text, rule_set.get("borderline", []))
                confidence = "borderline" if matched_rule else ""

            if not matched_rule:
                continue

            rows.append({
                "category": category,
                "confidence": confidence,
                "locus_tag": locus_tag,
                "old_locus_tag": old_locus_tag,
                "gene": gene,
                "product": product,
                "note": note,
                "contig": feature.get("contig", ""),
                "start": feature.get("start", ""),
                "end": feature.get("end", ""),
                "strand": feature.get("strand", ""),
                "matched_rule": matched_rule,
                "kleborate_support": kleb_reason,
            })

    result = pd.DataFrame(rows)

    if result.empty:
        return result

    # Deduplicate: same locus can match multiple categories in rare cases.
    # Keep all category assignments, but one row per category/locus_tag.
    confidence_rank = {"core": 0, "borderline": 1}
    result["confidence_rank"] = result["confidence"].map(confidence_rank).fillna(9)
    result = result.sort_values(["category", "locus_tag", "confidence_rank"])
    result = result.drop_duplicates(subset=["category", "locus_tag"], keep="first")
    result = result.drop(columns=["confidence_rank"])

    return result


# =============================================================================
# OUTPUT
# =============================================================================

def make_locus_tag_lists(extracted_df):
    if extracted_df.empty:
        return pd.DataFrame()

    lists = {}
    max_len = 0

    for category in sorted(extracted_df["category"].unique()):
        tags = (
            extracted_df[extracted_df["category"] == category]
            .sort_values(["confidence", "locus_tag"])["locus_tag"]
            .drop_duplicates()
            .tolist()
        )
        lists[category] = tags
        max_len = max(max_len, len(tags))

    padded = {}
    for category, tags in lists.items():
        padded[category] = tags + [""] * (max_len - len(tags))

    return pd.DataFrame(padded)


def make_summary(extracted_df):
    if extracted_df.empty:
        return pd.DataFrame(columns=["category", "confidence", "n_locus_tags"])

    return (
        extracted_df.groupby(["category", "confidence"])["locus_tag"]
        .nunique()
        .reset_index(name="n_locus_tags")
        .sort_values(["category", "confidence"])
    )


def write_output_excel(out_path, kleb_raw_df, kleb_modules_df, extracted_df, log_path):
    lists_df = make_locus_tag_lists(extracted_df)
    summary_df = make_summary(extracted_df)

    notes_df = pd.DataFrame([
        {
            "topic": "What Kleborate contributes",
            "note": (
                "Kleborate detects and types Klebsiella-specific modules such as ybt, clb, "
                "iuc, iro, rmp, K locus and O locus. The exact per-gene locus tags are then "
                "extracted from the supplied annotation file."
            ),
        },
        {
            "topic": "Current locus tags",
            "note": (
                "The script uses the current locus_tag field from the annotation file. "
                "old_locus_tag is reported separately when present but is not used as the main identifier."
            ),
        },
        {
            "topic": "Module-filtered families",
            "note": (
                "Yersiniabactin, aerobactin, salmochelin, colibactin and rmp categories are only "
                "extracted when the corresponding Kleborate module is detected, unless this behavior "
                "is disabled in USER SETTINGS."
            ),
        },
        {
            "topic": "Borderline hits",
            "note": (
                "Borderline hits are more general or context-dependent annotations. For strict "
                "manuscript analyses, inspect them manually."
            ),
        },
        {
            "topic": "Kleborate log",
            "note": str(log_path),
        },
    ])

    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        kleb_raw_df.to_excel(writer, sheet_name="Kleborate_raw", index=False)
        kleb_modules_df.to_excel(writer, sheet_name="Kleborate_modules", index=False)
        extracted_df.to_excel(writer, sheet_name="Extracted_genes", index=False)
        summary_df.to_excel(writer, sheet_name="Summary_counts", index=False)
        lists_df.to_excel(writer, sheet_name="Locus_tag_lists", index=False)
        notes_df.to_excel(writer, sheet_name="Notes", index=False)

        if not extracted_df.empty:
            for category in sorted(extracted_df["category"].unique()):
                sheet = safe_sheet_name(category)
                extracted_df[extracted_df["category"] == category].to_excel(
                    writer,
                    sheet_name=sheet,
                    index=False,
                )


# =============================================================================
# MAIN
# =============================================================================

def main():
    assembly_path, annotation_path, outdir = resolve_runtime_paths()
    outdir.mkdir(parents=True, exist_ok=True)

    sample_id = safe_sample_name(assembly_path)

    print(f"[INFO] Assembly:   {assembly_path}")
    print(f"[INFO] Annotation: {annotation_path}")
    print(f"[INFO] Output dir: {outdir}")

    annotation_df = parse_annotation(annotation_path)
    print(f"[INFO] Parsed annotation features: {len(annotation_df)}")

    if RUN_KLEBORATE:
        kleb_out, log_path = run_kleborate(assembly_path, outdir, sample_id)
        kleb_raw_df, kleb_table_path = read_kleborate_table(kleb_out)
    else:
        log_path = ""
        kleb_raw_df = pd.DataFrame()

    kleb_modules_df = detect_kleborate_modules(kleb_raw_df)
    print("[INFO] Kleborate modules:")
    print(kleb_modules_df)

    extracted_df = extract_genes_by_category(annotation_df, kleb_modules_df)

    if extracted_df.empty:
        print("[WARNING] No virulence-associated genes extracted.")
    else:
        print("[INFO] Extracted genes by category:")
        print(make_summary(extracted_df))

    out_path = outdir / OUTPUT_EXCEL_NAME
    write_output_excel(out_path, kleb_raw_df, kleb_modules_df, extracted_df, log_path)

    print(f"[DONE] Output written to:\n{out_path}")

    if messagebox is not None:
        try:
            messagebox.showinfo(
                "Done",
                f"Kleborate-guided Klebsiella virulence gene extraction completed.\n\n{out_path}"
            )
        except Exception:
            pass


if __name__ == "__main__":
    main()