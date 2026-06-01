#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ABRicate-based virulence functional screening pipeline

Purpose
-------
Screen assembled bacterial genomes for virulence-associated genes/loci using
ABRicate databases such as VFDB, VICTORS, Ecoli_VF and UPEC_ExPEC_VF.
Optional hooks are retained for VirulenceFinder and SPIFinder when these tools
are installed and configured locally.

The script merges parsed hits, assigns broad functional categories, optionally
maps hits to locus tags using matching annotation files, and exports summary
tables for downstream analysis.

Input metadata
--------------
A metadata file with at least these columns:

    sample_id,species,assembly

Optional annotation columns can be used for locus-tag mapping:

    annotation, annotation_file, gff, gff3, gbff, gbk, genbank

The metadata file can be CSV, TSV, semicolon-separated CSV, or Excel (.xlsx).

Minimal example:

    sample_id,species,assembly,annotation
    Ecoli_O157H7,ecoli_o157_h7,genomes/Ecoli_O157H7.fna,annotations/Ecoli_O157H7.gff
    SL1344,salmonella_typhimurium,genomes/SL1344.fna,annotations/SL1344.gff
    Sflex,shigella_flexneri,genomes/Shigella_flexneri.fna,annotations/Shigella_flexneri.gff

Outputs
-------
The main Excel workbook contains:

    - all_raw_hits
    - deduplicated_hits
    - functional_counts
    - functional_gene_lists
    - input_metadata
    - tool_run_log
    - expected_major_groups
    - notes

A clustering-ready workbook containing functional locus-tag lists can also be
generated when annotation files are provided.

Recommended use
---------------
Activate an environment containing the required command-line tools, then run:

    conda activate virulence_screen
    python GF11_abricate_virulence_functional_screen.py

If command-line arguments are omitted, the script opens file/directory pickers
to select the metadata file and output folder.
"""


# =============================================================================
# USER SETTINGS
# =============================================================================

from pathlib import Path

# ---------------------------------------------------------------------------
# Interactive / command-line input
# ---------------------------------------------------------------------------
# Leave these empty to use file/directory pickers at runtime.
# Alternatively, pass --metadata and --output-dir on the command line.
METADATA_FILE = ""
OUTPUT_DIR = ""

# Optional VirulenceFinder database path. Leave empty to select it interactively
# only when RUN_VIRULENCEFINDER=True.
VIRULENCEFINDER_DB_PATH = ""

# ---------------------------------------------------------------------------
# Tool switches
# ---------------------------------------------------------------------------
RUN_ABRICATE = True
RUN_VIRULENCEFINDER = False
RUN_SPIFINDER = False

# Optional SPIFinder command-line integration.
# Leave empty unless a working local SPIFinder command is available.
# Available placeholders: {assembly}, {outdir}, {sample}
SPIFINDER_CMD_TEMPLATE = ""

# ---------------------------------------------------------------------------
# ABRicate configuration
# ---------------------------------------------------------------------------
ABRICATE_BASE_DATABASES = ["vfdb", "victors"]
ABRICATE_ECOLI_SPECIFIC_DATABASES = ["ecoli_vf", "upec_expec_vf"]
ABRICATE_MIN_IDENTITY = 80
ABRICATE_MIN_COVERAGE = 60

# ---------------------------------------------------------------------------
# VirulenceFinder configuration
# ---------------------------------------------------------------------------
VIRULENCEFINDER_MIN_COVERAGE = 0.60
VIRULENCEFINDER_MIN_IDENTITY = 0.90

# ---------------------------------------------------------------------------
# General behavior
# ---------------------------------------------------------------------------
DRY_RUN_ONLY = False
CONTINUE_AFTER_TOOL_ERROR = True
WRITE_DEDUPLICATED_TABLE = True
SKIP_UNSUPPORTED_SPECIES = True

# Convert Windows paths such as C:\Users\... to WSL paths such as /mnt/c/Users/...
# when the script is launched from WSL/Linux.
AUTO_CONVERT_WINDOWS_PATHS_TO_WSL = True

CREATE_ABRICATE_SUMMARY = True

# ---------------------------------------------------------------------------
# Locus tag / Entrez Gene ID mapping
# ---------------------------------------------------------------------------
# ABRicate reports coordinates on contigs, but not current locus tags.
# To obtain correct locus tags, provide matching annotation files for the exact
# assemblies used here. Add an optional metadata column named annotation, gff,
# gff3, gbff, gbk, or genbank.
ENABLE_LOCUS_TAG_MAPPING = True
AUTO_DETECT_ANNOTATION_FILES = True
ANNOTATION_METADATA_COLUMNS = [
    "annotation", "annotation_file", "gff", "gff3", "gbff", "gbk", "genbank"
]
ANNOTATION_SUFFIXES = [".gff", ".gff3", ".gbff", ".gbk", ".gb", ".genbank"]
MIN_OVERLAP_FRACTION_OF_HIT = 0.50
MIN_OVERLAP_BP = 20

# Extra workbook for direct input to downstream clustering workflows.
CREATE_FUNCTIONAL_COUNTS_WORKBOOK = True
FUNCTIONAL_COUNTS_WORKBOOK_NAME = "Functional_counts.xlsx"


# =============================================================================
# IMPORTS
# =============================================================================

import os
import re
import io
import sys
import json
import shlex
import shutil
import argparse
import subprocess
from datetime import datetime

import pandas as pd

try:
    import tkinter as tk
    from tkinter import filedialog
except Exception:  # Tkinter may be unavailable on headless systems.
    tk = None
    filedialog = None

# Ensure command-line tools installed in the active conda environment are found,
# even when launched from Spyder using the virulence_screen interpreter.
CURRENT_ENV_BIN = str(Path(sys.executable).parent)
os.environ["PATH"] = CURRENT_ENV_BIN + os.pathsep + os.environ.get("PATH", "")


# =============================================================================
# SPECIES NORMALIZATION
# =============================================================================

SPECIES_ALIASES = {
    # EHEC / STEC O157:H7
    "e_coli_o157_h7": "ecoli_o157_h7",
    "escherichia_coli_o157_h7": "ecoli_o157_h7",
    "ecoli_o157_h7": "ecoli_o157_h7",
    "ecoli_o157": "ecoli_o157_h7",
    "e_coli_o157": "ecoli_o157_h7",
    "o157_h7": "ecoli_o157_h7",
    "o157": "ecoli_o157_h7",
    "ehec": "ecoli_o157_h7",
    "stec": "ecoli_o157_h7",
    "ehec_o157_h7": "ecoli_o157_h7",
    "stec_o157_h7": "ecoli_o157_h7",

    # UPEC UTI89
    "ecoli_uti89": "ecoli_upec_uti89",
    "e_coli_uti89": "ecoli_upec_uti89",
    "escherichia_coli_uti89": "ecoli_upec_uti89",
    "upec_uti89": "ecoli_upec_uti89",
    "uti89": "ecoli_upec_uti89",
    "upec": "ecoli_upec",

    # Generic E. coli, if needed later
    "ecoli": "ecoli_generic",
    "e_coli": "ecoli_generic",
    "escherichia_coli": "ecoli_generic",
    "escherichia": "ecoli_generic",

    # Salmonella
    "salmonella": "salmonella_enterica",
    "salmonella_enterica": "salmonella_enterica",
    "salmonella_typhimurium": "salmonella_typhimurium",
    "salmonella_enterica_serovar_typhimurium": "salmonella_typhimurium",
    "s_typhimurium": "salmonella_typhimurium",
    "sl1344": "salmonella_typhimurium",

    # Shigella
    "shigella": "shigella_flexneri",
    "shigella_flexneri": "shigella_flexneri",
    "shigella_flex": "shigella_flexneri",
    "s_flexneri": "shigella_flexneri",
    "sflex": "shigella_flexneri",

    # Klebsiella
    "klebsiella": "klebsiella_pneumoniae",
    "klebsiella_pneumoniae": "klebsiella_pneumoniae",
    "klebsiella_pneum": "klebsiella_pneumoniae",
    "k_pneumoniae": "klebsiella_pneumoniae",
    "kp": "klebsiella_pneumoniae",

    # Enterobacter cloacae complex
    "enterobacter": "enterobacter_cloacae",
    "enterobacter_cloacae": "enterobacter_cloacae",
    "enterobacter_cloacae_complex": "enterobacter_cloacae",
    "e_cloacae": "enterobacter_cloacae",
    "enterobacter_cloaceae": "enterobacter_cloacae",  # common typo
    "ecloaceae": "enterobacter_cloacae",
    "ecloa": "enterobacter_cloacae",
    "eclo": "enterobacter_cloacae",

    # Unsupported / intentionally removed from this Enterobacteriaceae-focused analysis
    "enterococcus_faecium": "unsupported_enterococcus_faecium",
    "e_faecium": "unsupported_enterococcus_faecium",
}

ECOLI_LIKE_SPECIES = {
    "ecoli_o157_h7",
    "ecoli_upec_uti89",
    "ecoli_upec",
    "ecoli_generic",
}

SHIGELLA_LIKE_SPECIES = {
    "shigella_flexneri",
}

SALMONELLA_LIKE_SPECIES = {
    "salmonella_enterica",
    "salmonella_typhimurium",
}

KLEBSIELLA_LIKE_SPECIES = {
    "klebsiella_pneumoniae",
}

ENTEROBACTER_LIKE_SPECIES = {
    "enterobacter_cloacae",
}

SUPPORTED_SPECIES = (
    ECOLI_LIKE_SPECIES
    | SHIGELLA_LIKE_SPECIES
    | SALMONELLA_LIKE_SPECIES
    | KLEBSIELLA_LIKE_SPECIES
    | ENTEROBACTER_LIKE_SPECIES
)


def normalize_species_name(value: str) -> str:
    """Normalize user-provided species/pathotype names to internal identifiers."""
    original = str(value).strip()
    s = original.lower()
    s = s.replace(":", " ")
    s = re.sub(r"[^a-z0-9]+", "_", s).strip("_")

    if s in SPECIES_ALIASES:
        return SPECIES_ALIASES[s]

    # Fuzzy fallbacks for truncated spreadsheet labels or informal names.
    if "uti89" in s or ("upec" in s and "coli" in s):
        return "ecoli_upec_uti89"
    if "o157" in s or "ehec" in s or "stec" in s:
        return "ecoli_o157_h7"
    if "salmonella" in s or "typhimurium" in s or "sl1344" in s:
        return "salmonella_typhimurium" if ("typhimurium" in s or "sl1344" in s) else "salmonella_enterica"
    if "shigella" in s or "flexneri" in s or "sflex" in s:
        return "shigella_flexneri"
    if "klebsiella" in s or "pneumoniae" in s:
        return "klebsiella_pneumoniae"
    if "enterobacter" in s or "cloaca" in s or "cloaceae" in s or "ecloa" in s:
        return "enterobacter_cloacae"

    return s


def is_ecoli_or_shigella_like(species: str) -> bool:
    return species in ECOLI_LIKE_SPECIES or species in SHIGELLA_LIKE_SPECIES


def is_salmonella_like(species: str) -> bool:
    return species in SALMONELLA_LIKE_SPECIES


# =============================================================================
# FUNCTIONAL GROUP CLASSIFIER
# =============================================================================

def _contains_any(text: str, patterns) -> bool:
    return any(re.search(p, text, flags=re.IGNORECASE) for p in patterns)


def classify_virulence_hit(gene_name: str, product: str = "", source_text: str = "") -> str:
    """
    Rule-based functional grouping.

    This is intentionally broad. It is a biologically interpretable first-pass
    classification, not a substitute for manual curation before publication.
    """
    gene = str(gene_name or "").strip()
    product = str(product or "").strip()
    text = f"{gene} {product} {source_text}".lower()

    # --- Toxins and genotoxins ---
    if _contains_any(text, [
        r"\bstx[12]?[ab]?\b", r"shiga", r"verotoxin", r"vero toxin",
        r"\bsubab\b", r"shiga-like",
    ]):
        return "Toxin - Shiga/verotoxin"

    if _contains_any(text, [
        r"\bclb[a-q]\b", r"colibactin", r"genotoxin", r"\bcdt[abc]\b",
        r"cytolethal distending toxin",
    ]):
        return "Toxin - genotoxin"

    if _contains_any(text, [
        r"\bhly[a-z]?\b", r"hemolysin", r"haemolysin",
        r"\bcyl[a-z]?\b", r"cytolysin", r"\bcnf[0-9]?\b",
        r"\bastA\b", r"\best[a-z0-9]?\b", r"\belt[a-z0-9]?\b",
        r"enterotoxin", r"heat-labile toxin", r"heat-stable toxin",
    ]):
        return "Toxin - other"

    # --- T3SS effectors ---
    if _contains_any(text, [
        # EHEC/EPEC LEE and non-LEE effectors
        r"\btir\b", r"\bmap\b", r"\bcif\b", r"\bnle[a-z0-9_]*\b",
        r"\besp[fghijklmnpqrstuvwxz][0-9a-z_]*\b", r"\btccp\b",

        # Salmonella effectors
        r"\bsop[a-z0-9_]*\b", r"\bsptp\b", r"\bavra\b",
        r"\bsif[a-z0-9_]*\b", r"\bsse[a-z0-9_]*\b",
        r"\bpip[a-z0-9_]*\b", r"\bste[a-z0-9_]*\b",
        r"\bgog[a-z0-9_]*\b", r"\bgtg[a-z0-9_]*\b",
        r"\bslrp\b", r"\bsar[a-z0-9_]*\b",

        # Shigella effectors
        r"\bipa[h][0-9a-z_]*\b", r"\bosp[a-z0-9_]*\b",
        r"\bvir[a-z0-9_]*\b.*effector",
    ]):
        return "Type III secretion - effector"

    # --- T3SS apparatus / translocon / chaperones ---
    if _contains_any(text, [
        r"type iii", r"type 3 secretion", r"t3ss", r"injectisome", r"\bsct[a-z]\b",

        # EHEC/EPEC LEE apparatus/translocon
        r"\besc[a-z]\b", r"\bsep[a-z]\b", r"\bces[a-z]\b",
        r"\besp[abd]\b", r"\beae\b", r"intimin",

        # Salmonella SPI-1 apparatus/translocon/chaperones
        r"\binv[a-z]\b", r"\bprg[a-z]\b", r"\borg[a-z]\b",
        r"\bsip[abcd]\b", r"\bsic[ap]\b", r"\bspa[opqrs]\b",

        # Salmonella SPI-2 apparatus/chaperones
        r"\bssa[a-z]\b", r"\bssc[a-z]\b", r"\bssr[ab]\b",

        # Shigella Mxi-Spa apparatus / Ipa translocon
        r"\bmxi[a-z]\b", r"\bspa[0-9a-z]?\b", r"\bipa[abcd]\b", r"\bipg[a-z]\b",
    ]):
        return "Type III secretion - apparatus/translocon"

    # --- T6SS ---
    if _contains_any(text, [
        r"type vi", r"type 6 secretion", r"t6ss",
        r"\btss[a-z]\b", r"\bhcp\b", r"\bvgrg\b", r"\bclpv\b",
        r"\bvas[a-z]\b", r"\bimp[a-z]\b", r"rhs toxin",
    ]):
        return "Type VI secretion / interbacterial competition"

    # --- Siderophores and iron ---
    if _contains_any(text, [
        r"\bybt[a-z]\b", r"\bfyua\b", r"\birp[12]\b", r"yersiniabactin",
        r"\biuc[a-d]\b", r"\biuta\b", r"aerobactin",
        r"\biro[bcden]\b", r"salmochelin",
        r"\bent[abcdefs]\b", r"\bfep[a-g]\b", r"\bfes\b", r"enterobactin",
        r"\bfec[a-e]\b", r"\bfhu[a-f]\b", r"\bsit[abcd]\b", r"\bsit\b",
        r"\bchu[a-z]\b", r"\btonb\b", r"iron acquisition", r"siderophore",
        r"heme uptake", r"ferric", r"ferrous",
    ]):
        return "Iron acquisition / siderophore"

    # --- Capsule and surface polysaccharide ---
    if _contains_any(text, [
        r"capsule", r"capsular", r"\bcps\b", r"\bk_locus\b", r"\bk locus\b",
        r"\bwzi\b", r"\bwza\b", r"\bwzb\b", r"\bwzc\b",
        r"\bgalf\b", r"\bugd\b", r"\bwca[a-z]\b", r"colanic acid",
        r"exopolysaccharide", r"\brmpa2?\b", r"\brmpd\b", r"hypermucovisc",
        r"hypermucoviscosity", r"mucoviscosity",
    ]):
        return "Capsule / surface polysaccharide"

    # --- LPS and O antigen ---
    if _contains_any(text, [
        r"lipopolysaccharide", r"\blps\b", r"o-antigen", r"\bo_locus\b", r"\bo locus\b",
        r"\bwaa[a-z]\b", r"\brfa[a-z]\b", r"\brfb[a-z]\b",
        r"\bwzx\b", r"\bwzy\b", r"\brml[a-d]\b",
    ]):
        return "LPS / O-antigen"

    # --- Adhesins, pili, fimbriae ---
    if _contains_any(text, [
        r"adhesin", r"adhesion", r"fimbr", r"\bfim[a-z]\b", r"\bmrk[a-z]\b",
        r"\bpap[a-z]\b", r"\bsfa[a-z]\b", r"\bfoc[a-z]\b",
        r"\bafa[a-z]?\b", r"\bdaa[a-z]?\b", r"\byad[a-z]\b",
        r"pilus", r"pili", r"\bpil[a-z]\b",
        r"\bcsg[a-g]\b", r"curli",
        r"\beae\b", r"intimin", r"\biha\b", r"\bsaa\b", r"\btoxB\b",
        r"autotransporter", r"\bag43\b", r"\bflu\b",
    ]):
        return "Adhesin / fimbriae / pilus"

    # --- Biofilm ---
    if _contains_any(text, [
        r"biofilm", r"\bpga[abcd]\b", r"\bbcs[a-z]\b", r"\bbap\b", r"cellulose",
    ]):
        return "Biofilm-associated"

    # --- Motility / cell-to-cell spread ---
    if _contains_any(text, [
        r"\bicsa\b", r"\bvirg\b", r"actin.*motility",
        r"cell-to-cell spread", r"intercellular spread",
        r"\bfli[cgdm]\b", r"\bflg[a-z]\b", r"\bflh[a-z]\b",
        r"flagell", r"\bmot[ab]\b", r"swarming",
    ]):
        return "Motility / cell-to-cell spread"

    # --- Intracellular survival / stress ---
    if _contains_any(text, [
        r"\bmgtc\b", r"\bpho[pq]\b", r"\bpag[a-z]\b", r"\bmig[a-z]\b",
        r"\bspv[abcd]\b", r"\bsod[abc]\b", r"\bkat[eg]\b",
        r"acid resistance", r"\bgad[abcex]\b", r"\badi[ac]\b",
        r"intracellular survival", r"macrophage survival", r"oxidative stress",
        r"nitrosative stress",
    ]):
        return "Intracellular survival / stress resistance"

    # --- Immune evasion / serum resistance ---
    if _contains_any(text, [
        r"serum resistance", r"\biss\b", r"\btrat\b", r"\bompA\b",
        r"immune evasion", r"complement resistance", r"factor h",
    ]):
        return "Immune evasion / serum resistance"

    # --- Virulence regulation ---
    if _contains_any(text, [
        r"\bler\b", r"\bgrl[ar]\b", r"\bhil[acd]\b", r"\binv[f]\b",
        r"\bssr[ab]\b", r"\bpho[pq]\b", r"\bslyA\b",
        r"\bvir[fbr]\b", r"\brmpa2?\b", r"virulence regulator",
        r"transcriptional regulator.*virulence",
    ]):
        return "Virulence regulation"

    # --- Plasmid / pathogenicity island generic ---
    if _contains_any(text, [
        r"pathogenicity island", r"\bspi[-_ ]?[0-9]+\b", r"\bpai\b",
        r"virulence plasmid", r"pseudogene.*virulence",
    ]):
        return "Pathogenicity island / virulence plasmid"

    return "Unknown virulence-associated"


# =============================================================================
# UTILITIES
# =============================================================================

def now_str():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def safe_sample_id(value: str) -> str:
    s = str(value).strip()
    s = re.sub(r"[^A-Za-z0-9_.-]+", "_", s)
    return s.strip("_")


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)


def tool_exists(executable: str) -> bool:
    return shutil.which(executable) is not None


def _tk_available() -> bool:
    return tk is not None and filedialog is not None


def choose_file(title: str, patterns):
    """Open a file picker and return the selected path."""
    if not _tk_available():
        raise RuntimeError(
            "Tkinter file picker is unavailable. Provide paths with command-line "
            "arguments, for example: --metadata input_metadata.csv --output-dir results"
        )

    root = tk.Tk()
    root.withdraw()
    path = filedialog.askopenfilename(title=title, filetypes=patterns)
    root.destroy()

    if not path:
        raise SystemExit(f"No file selected for: {title}")
    return Path(path)


def choose_directory(title: str):
    """Open a directory picker and return the selected path."""
    if not _tk_available():
        raise RuntimeError(
            "Tkinter directory picker is unavailable. Provide --output-dir on the command line."
        )

    root = tk.Tk()
    root.withdraw()
    path = filedialog.askdirectory(title=title)
    root.destroy()

    if not path:
        raise SystemExit(f"No directory selected for: {title}")
    return Path(path)


def parse_command_line_args():
    parser = argparse.ArgumentParser(
        description="ABRicate-based virulence functional screen for assembled bacterial genomes."
    )
    parser.add_argument(
        "--metadata",
        default=METADATA_FILE,
        help="Input metadata file with columns sample_id, species and assembly.",
    )
    parser.add_argument(
        "--output-dir",
        default=OUTPUT_DIR,
        help="Output directory for Excel/TSV results and raw tool outputs.",
    )
    parser.add_argument(
        "--virulencefinder-db",
        default=VIRULENCEFINDER_DB_PATH,
        help="Optional VirulenceFinder database directory.",
    )
    parser.add_argument(
        "--no-gui",
        action="store_true",
        help="Disable file pickers and require paths to be supplied as arguments or USER SETTINGS.",
    )
    return parser.parse_args()


def resolve_runtime_paths():
    """Resolve metadata, output and optional database paths from CLI, settings or pickers."""
    args = parse_command_line_args()

    metadata_file = str(args.metadata).strip()
    output_dir = str(args.output_dir).strip()
    vf_db_path = str(args.virulencefinder_db).strip()

    if not metadata_file:
        if args.no_gui:
            raise SystemExit("Missing --metadata. Provide a metadata CSV/TSV/XLSX file.")
        metadata_file = str(choose_file(
            "Select metadata file",
            [
                ("Metadata files", "*.csv *.tsv *.txt *.xlsx *.xlsm *.xls"),
                ("All files", "*.*"),
            ],
        ))

    if not output_dir:
        if args.no_gui:
            raise SystemExit("Missing --output-dir. Provide an output directory.")
        output_dir = str(choose_directory("Select output directory"))

    return Path(metadata_file), Path(output_dir), str(vf_db_path)


def windows_path_to_wsl(path_text: str) -> str:
    """Convert C:\\something to /mnt/c/something when running under WSL/Linux."""
    p = str(path_text).strip().strip('"').strip("'")
    if not AUTO_CONVERT_WINDOWS_PATHS_TO_WSL:
        return p

    # Only convert on non-Windows platforms.
    if os.name == "nt":
        return p

    m = re.match(r"^([A-Za-z]):[\\/](.*)$", p)
    if not m:
        return p

    drive = m.group(1).lower()
    rest = m.group(2).replace("\\", "/")
    return f"/mnt/{drive}/{rest}"


def run_command(cmd, stdout_path=None, log_path=None, env=None, shell=False):
    """Run a command and return a structured log dictionary."""
    start = now_str()
    cmd_display = cmd if isinstance(cmd, str) else " ".join(shlex.quote(str(x)) for x in cmd)

    log_record = {
        "time": start,
        "command": cmd_display,
        "returncode": None,
        "stdout_file": str(stdout_path) if stdout_path else "",
        "log_file": str(log_path) if log_path else "",
        "status": "not_run",
        "message": "",
    }

    if DRY_RUN_ONLY:
        log_record["status"] = "dry_run"
        log_record["message"] = "Command not executed because DRY_RUN_ONLY=True."
        return log_record

    try:
        result = subprocess.run(
            cmd,
            shell=shell,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env,
        )

        log_record["returncode"] = result.returncode
        log_record["status"] = "ok" if result.returncode == 0 else "error"

        if stdout_path:
            Path(stdout_path).write_text(result.stdout, encoding="utf-8", errors="replace")

        if log_path:
            with open(log_path, "w", encoding="utf-8", errors="replace") as f:
                f.write(f"[TIME]\n{start}\n\n")
                f.write(f"[COMMAND]\n{cmd_display}\n\n")
                f.write(f"[RETURNCODE]\n{result.returncode}\n\n")
                f.write("[STDOUT]\n")
                f.write(result.stdout or "")
                f.write("\n\n[STDERR]\n")
                f.write(result.stderr or "")

        if result.returncode != 0:
            log_record["message"] = (result.stderr or result.stdout or "").strip()[:2000]
            if not CONTINUE_AFTER_TOOL_ERROR:
                raise RuntimeError(f"Command failed: {cmd_display}\n{log_record['message']}")

        return log_record

    except Exception as exc:
        log_record["status"] = "exception"
        log_record["message"] = str(exc)
        if log_path:
            Path(log_path).write_text(
                f"[TIME]\n{start}\n\n[COMMAND]\n{cmd_display}\n\n[EXCEPTION]\n{exc}\n",
                encoding="utf-8",
                errors="replace",
            )
        if not CONTINUE_AFTER_TOOL_ERROR:
            raise
        return log_record


def read_table_from_first_header(path: Path):
    """Read a TSV/CSV-like file even if the tool added preamble lines."""
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()

    text = path.read_text(encoding="utf-8", errors="replace")
    lines = [ln for ln in text.splitlines() if ln.strip()]

    if not lines:
        return pd.DataFrame()

    header_idx = None
    header_markers = [
        "#FILE", "FILE\t", "SEQUENCE\t", "Gene\t", "GENE\t",
        "gene\t", "Sample\t", "sample\t", "Assembly\t", "assembly\t",
        "strain\t", "species\t", "input_file\t"
    ]

    for i, ln in enumerate(lines):
        if any(ln.startswith(m) for m in header_markers):
            header_idx = i
            break

    content = "\n".join(lines[header_idx:]) if header_idx is not None else "\n".join(lines)
    first_line = content.splitlines()[0]

    if "\t" in first_line:
        sep = "\t"
    elif ";" in first_line and first_line.count(";") > first_line.count(","):
        sep = ";"
    else:
        sep = ","

    try:
        df = pd.read_csv(io.StringIO(content), sep=sep, dtype=str)
    except Exception:
        return pd.DataFrame()

    df.columns = [str(c).strip().lstrip("#") for c in df.columns]
    return df.fillna("")


def to_float_or_blank(value):
    try:
        if value is None or str(value).strip() == "":
            return ""
        return float(str(value).replace("%", "").strip())
    except Exception:
        return ""


def standard_hit_row(
    sample_id,
    species,
    tool,
    database,
    gene_name="",
    product="",
    contig="",
    start="",
    end="",
    strand="",
    identity="",
    coverage="",
    accession="",
    evidence="homology_or_tool_call",
    source_file="",
    extra="",
):
    gene_name = str(gene_name or "").strip()
    product = str(product or "").strip()
    extra = str(extra or "").strip()

    functional_group = classify_virulence_hit(gene_name, product, extra)

    return {
        "sample_id": sample_id,
        "species": species,
        "tool": tool,
        "database": database,
        "gene_name": gene_name,
        "product_or_locus": product,
        "functional_group": functional_group,
        "contig": contig,
        "start": start,
        "end": end,
        "strand": strand,
        "identity_percent": to_float_or_blank(identity),
        "coverage_percent": to_float_or_blank(coverage),
        "accession": accession,
        "evidence": evidence,
        "source_file": str(source_file),
        "extra": extra,
    }


# =============================================================================
# ABRICATE DATABASE HELPERS
# =============================================================================

def get_available_abricate_databases():
    """Return a set of ABRicate database names available in the current environment."""
    if not tool_exists("abricate"):
        return set()

    try:
        result = subprocess.run(
            ["abricate", "--list"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
    except Exception:
        return set()

    if result.returncode != 0:
        return set()

    dbs = set()
    for line in result.stdout.splitlines():
        if not line.strip() or line.startswith("DATABASE"):
            continue
        dbs.add(line.split()[0])
    return dbs


# =============================================================================
# ANNOTATION / LOCUS TAG MAPPING
# =============================================================================

def normalize_contig_id(value: str) -> str:
    """Normalize contig/record IDs for matching ABRicate hits to annotations."""
    s = str(value or "").strip()
    if not s:
        return ""
    s = s.lstrip(">")
    # ABRicate and GFF/GenBank usually use the first token before whitespace.
    s = s.split()[0]
    return s


def parse_dbxref_for_entrez(value: str) -> str:
    """Extract Entrez GeneID from GFF/GenBank Dbxref/db_xref strings."""
    text = str(value or "")
    matches = re.findall(r"(?:GeneID|geneid|EntrezGene|Entrez):?\s*([0-9]+)", text)
    if matches:
        return ";".join(sorted(set(matches), key=matches.index))
    return ""


def parse_gff_attributes(attr_text: str) -> dict:
    """Parse a GFF3/GTF-like attributes column into a dictionary."""
    attrs = {}
    for part in str(attr_text or "").strip().split(";"):
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
        attrs[key.strip()] = value.strip()
    return attrs


def _first_nonempty(*values):
    for value in values:
        if value is None:
            continue
        if isinstance(value, (list, tuple)):
            if value:
                value = value[0]
            else:
                continue
        value = str(value).strip()
        if value:
            return value
    return ""


def _feature_record(contig, feature_type, start, end, strand, locus_tag, entrez_gene_id,
                    gene_name, product, source, annotation_file):
    return {
        "contig": normalize_contig_id(contig),
        "feature_type": str(feature_type or ""),
        "feature_start": int(start),
        "feature_end": int(end),
        "feature_strand": str(strand or ""),
        "locus_tag": str(locus_tag or ""),
        "EntrezGeneID": str(entrez_gene_id or ""),
        "mapped_gene_name": str(gene_name or ""),
        "mapped_product": str(product or ""),
        "annotation_source": str(source or ""),
        "annotation_file": str(annotation_file or ""),
    }


def parse_gff_annotation(annotation_file: Path):
    """Parse CDS/gene features from a GFF/GFF3 file."""
    features = []
    annotation_file = Path(annotation_file)

    with open(annotation_file, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue

            seqid, source, feature_type, start, end, score, strand, phase, attrs_text = parts[:9]
            feature_type_l = feature_type.lower()
            if feature_type_l not in {"cds", "gene", "pseudogene"}:
                continue

            attrs = parse_gff_attributes(attrs_text)
            dbxref = _first_nonempty(
                attrs.get("Dbxref"), attrs.get("db_xref"), attrs.get("dbxref"), attrs.get("Note")
            )
            entrez_gene_id = parse_dbxref_for_entrez(dbxref)

            locus_tag = _first_nonempty(
                attrs.get("locus_tag"),
                attrs.get("old_locus_tag") if not attrs.get("locus_tag") else "",
                attrs.get("ID") if feature_type_l == "cds" else "",
            )

            gene_name = _first_nonempty(attrs.get("gene"), attrs.get("Name"), attrs.get("gene_name"))
            product = _first_nonempty(attrs.get("product"), attrs.get("Note"), attrs.get("function"))

            try:
                features.append(_feature_record(
                    seqid, feature_type, int(start), int(end), strand, locus_tag, entrez_gene_id,
                    gene_name, product, source, annotation_file
                ))
            except Exception:
                continue

    return features


def parse_genbank_annotation(annotation_file: Path):
    """Parse CDS/gene features from a GenBank/GBFF/GBK file using Biopython."""
    try:
        from Bio import SeqIO
    except Exception as exc:
        print(f"[WARNING] Biopython is required to parse GenBank files but could not be imported: {exc}")
        return []

    features = []
    annotation_file = Path(annotation_file)

    try:
        records = SeqIO.parse(str(annotation_file), "genbank")
        for record in records:
            contig = normalize_contig_id(record.id or record.name)
            for feat in record.features:
                if feat.type.lower() not in {"cds", "gene", "pseudogene"}:
                    continue
                quals = feat.qualifiers
                locus_tag = _first_nonempty(quals.get("locus_tag"), quals.get("old_locus_tag"))
                gene_name = _first_nonempty(quals.get("gene"), quals.get("gene_synonym"))
                product = _first_nonempty(quals.get("product"), quals.get("note"), quals.get("function"))
                dbxref_text = ";".join(quals.get("db_xref", []))
                entrez_gene_id = parse_dbxref_for_entrez(dbxref_text)
                strand = "+" if feat.location.strand == 1 else "-" if feat.location.strand == -1 else ""
                start = int(feat.location.start) + 1
                end = int(feat.location.end)

                features.append(_feature_record(
                    contig, feat.type, start, end, strand, locus_tag, entrez_gene_id,
                    gene_name, product, "GenBank", annotation_file
                ))
    except Exception as exc:
        print(f"[WARNING] Could not parse GenBank annotation file {annotation_file}: {exc}")
        return []

    return features


def parse_annotation_file(annotation_file):
    """Parse GFF/GFF3 or GenBank annotation and return a list of feature records."""
    if not annotation_file:
        return []
    annotation_file = Path(annotation_file)
    if not annotation_file.exists():
        return []

    suffix = annotation_file.suffix.lower()
    if suffix in {".gff", ".gff3"}:
        return parse_gff_annotation(annotation_file)
    if suffix in {".gbff", ".gbk", ".gb", ".genbank"}:
        return parse_genbank_annotation(annotation_file)

    # Try GFF-like parsing first, then GenBank as fallback.
    gff_features = parse_gff_annotation(annotation_file)
    if gff_features:
        return gff_features
    return parse_genbank_annotation(annotation_file)


def build_annotation_index(features):
    """Index features by normalized contig ID."""
    index = {}
    for feat in features:
        contig = normalize_contig_id(feat.get("contig", ""))
        if not contig:
            continue
        index.setdefault(contig, []).append(feat)

    for contig in index:
        index[contig].sort(key=lambda x: (x["feature_start"], x["feature_end"]))

    return index


def find_best_feature_for_hit(contig, start, end, annotation_index):
    """Find the best annotation feature overlapping an ABRicate-like hit."""
    contig_norm = normalize_contig_id(contig)
    if not contig_norm:
        return None, "missing_contig", {}

    if contig_norm not in annotation_index:
        return None, "contig_not_found_in_annotation", {"normalized_contig": contig_norm}

    try:
        hit_start = int(float(start))
        hit_end = int(float(end))
    except Exception:
        return None, "missing_or_invalid_coordinates", {"normalized_contig": contig_norm}

    if hit_start > hit_end:
        hit_start, hit_end = hit_end, hit_start

    hit_len = max(1, hit_end - hit_start + 1)
    candidates = []

    for feat in annotation_index.get(contig_norm, []):
        f_start = int(feat["feature_start"])
        f_end = int(feat["feature_end"])
        overlap = max(0, min(hit_end, f_end) - max(hit_start, f_start) + 1)
        if overlap <= 0:
            continue

        frac_hit = overlap / hit_len
        frac_feature = overlap / max(1, f_end - f_start + 1)
        if overlap < MIN_OVERLAP_BP or frac_hit < MIN_OVERLAP_FRACTION_OF_HIT:
            continue

        is_cds = 1 if str(feat.get("feature_type", "")).lower() == "cds" else 0
        has_locus = 1 if str(feat.get("locus_tag", "")).strip() else 0
        score = (has_locus, is_cds, frac_hit, overlap, frac_feature)
        candidates.append((score, feat, {
            "overlap_bp": overlap,
            "overlap_fraction_of_hit": frac_hit,
            "overlap_fraction_of_feature": frac_feature,
            "normalized_contig": contig_norm,
        }))

    if not candidates:
        return None, "no_overlapping_feature_passing_threshold", {"normalized_contig": contig_norm}

    candidates.sort(key=lambda x: x[0], reverse=True)
    return candidates[0][1], "coordinate_overlap", candidates[0][2]


def annotate_hits_with_locus_tags(hits_df: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    """Add locus_tag and EntrezGeneID columns to hit tables using annotation files."""
    if hits_df.empty:
        return hits_df

    df = hits_df.copy()

    new_cols = {
        "locus_tag": "",
        "EntrezGeneID": "",
        "locus_mapping_status": "not_attempted",
        "locus_mapping_method": "",
        "annotation_file": "",
        "mapped_feature_type": "",
        "mapped_gene_name": "",
        "mapped_product": "",
        "mapped_feature_start": "",
        "mapped_feature_end": "",
        "mapped_feature_strand": "",
        "overlap_bp": "",
        "overlap_fraction_of_hit": "",
        "overlap_fraction_of_feature": "",
    }
    for col, default in new_cols.items():
        if col not in df.columns:
            df[col] = default

    if not ENABLE_LOCUS_TAG_MAPPING:
        df["locus_mapping_status"] = "disabled"
        return df

    metadata_by_sample = metadata.set_index("sample_id", drop=False).to_dict(orient="index")
    annotation_cache = {}

    for idx, row in df.iterrows():
        sample_id = row.get("sample_id", "")
        meta = metadata_by_sample.get(sample_id)
        if meta is None:
            df.at[idx, "locus_mapping_status"] = "sample_not_found_in_metadata"
            continue

        annotation_file = str(meta.get("annotation_file", "") or "").strip()
        if not annotation_file:
            df.at[idx, "locus_mapping_status"] = "no_annotation_file"
            continue

        annotation_path = Path(annotation_file)
        if not annotation_path.exists():
            df.at[idx, "locus_mapping_status"] = "annotation_file_not_found"
            df.at[idx, "annotation_file"] = str(annotation_path)
            continue

        df.at[idx, "annotation_file"] = str(annotation_path)

        if annotation_path not in annotation_cache:
            features = parse_annotation_file(annotation_path)
            annotation_cache[annotation_path] = build_annotation_index(features)

        annotation_index = annotation_cache[annotation_path]
        if not annotation_index:
            df.at[idx, "locus_mapping_status"] = "annotation_parsed_no_features"
            continue

        feature, status, overlap_info = find_best_feature_for_hit(
            row.get("contig", ""), row.get("start", ""), row.get("end", ""), annotation_index
        )

        df.at[idx, "locus_mapping_status"] = status
        df.at[idx, "locus_mapping_method"] = "coordinate_overlap" if feature else ""

        if feature:
            df.at[idx, "locus_tag"] = feature.get("locus_tag", "")
            df.at[idx, "EntrezGeneID"] = feature.get("EntrezGeneID", "")
            df.at[idx, "mapped_feature_type"] = feature.get("feature_type", "")
            df.at[idx, "mapped_gene_name"] = feature.get("mapped_gene_name", "")
            df.at[idx, "mapped_product"] = feature.get("mapped_product", "")
            df.at[idx, "mapped_feature_start"] = feature.get("feature_start", "")
            df.at[idx, "mapped_feature_end"] = feature.get("feature_end", "")
            df.at[idx, "mapped_feature_strand"] = feature.get("feature_strand", "")
            df.at[idx, "overlap_bp"] = overlap_info.get("overlap_bp", "")
            df.at[idx, "overlap_fraction_of_hit"] = round(float(overlap_info.get("overlap_fraction_of_hit", 0)), 4)
            df.at[idx, "overlap_fraction_of_feature"] = round(float(overlap_info.get("overlap_fraction_of_feature", 0)), 4)
        else:
            if "normalized_contig" in overlap_info:
                df.at[idx, "locus_mapping_method"] = f"normalized_contig={overlap_info['normalized_contig']}"

    return df


def find_annotation_file_for_assembly(assembly_path: str, row: pd.Series = None) -> str:
    """Find an annotation file from metadata columns or by basename autodetection."""
    if row is not None:
        lower_cols = {str(c).lower().strip(): c for c in row.index}
        for candidate in ANNOTATION_METADATA_COLUMNS:
            if candidate in lower_cols:
                value = str(row.get(lower_cols[candidate], "") or "").strip()
                if value:
                    value = windows_path_to_wsl(value)
                    if Path(value).exists():
                        return str(Path(value))
                    return str(Path(value))

    if not AUTO_DETECT_ANNOTATION_FILES:
        return ""

    assembly = Path(assembly_path)
    folder = assembly.parent
    stem = assembly.stem
    # Handles .fna.gz/.fa.gz style names by removing the second suffix.
    if assembly.suffix.lower() == ".gz":
        stem = Path(stem).stem

    candidates = []
    for suffix in ANNOTATION_SUFFIXES:
        candidates.append(folder / f"{stem}{suffix}")

    # Also check a common annotations subfolder.
    for suffix in ANNOTATION_SUFFIXES:
        candidates.append(folder / "annotations" / f"{stem}{suffix}")
        candidates.append(folder.parent / "annotations" / f"{stem}{suffix}")

    for candidate in candidates:
        if candidate.exists():
            return str(candidate)

    return ""


# =============================================================================
# PARSERS
# =============================================================================

def parse_abricate_output(path, sample_id, species, database):
    df = read_table_from_first_header(Path(path))
    rows = []

    if df.empty:
        return rows

    colmap = {c.upper(): c for c in df.columns}

    def get(row, col, default=""):
        c = colmap.get(col.upper())
        return row.get(c, default) if c else default

    for _, row in df.iterrows():
        gene = get(row, "GENE")
        product = get(row, "PRODUCT")
        if not gene and not product:
            continue

        rows.append(standard_hit_row(
            sample_id=sample_id,
            species=species,
            tool="ABRicate",
            database=database,
            gene_name=gene,
            product=product,
            contig=get(row, "SEQUENCE"),
            start=get(row, "START"),
            end=get(row, "END"),
            strand=get(row, "STRAND"),
            identity=get(row, "%IDENTITY"),
            coverage=get(row, "%COVERAGE"),
            accession=get(row, "ACCESSION"),
            evidence="ABRicate nucleotide homology hit",
            source_file=path,
            extra=(
                f"coverage_map={get(row, 'COVERAGE_MAP')}; "
                f"gaps={get(row, 'GAPS')}; "
                f"resistance={get(row, 'RESISTANCE')}"
            ),
        ))

    return rows


def parse_generic_table_outputs(root_dir, sample_id, species, tool_name, database_name=""):
    """Flexible parser for TSV/CSV outputs produced by optional tools."""
    root_dir = Path(root_dir)
    rows = []

    if not root_dir.exists():
        return rows

    if root_dir.is_file():
        files = [root_dir]
    else:
        allowed_suffixes = {".tsv", ".tab", ".txt", ".csv"}
        files = [p for p in root_dir.rglob("*") if p.is_file() and p.suffix.lower() in allowed_suffixes]

    for p in files:
        if p.stat().st_size == 0 or p.stat().st_size > 50_000_000:
            continue
        if "log" in p.name.lower():
            continue

        df = read_table_from_first_header(p)
        if df.empty or len(df.columns) < 2:
            continue

        lower_cols = {c.lower().strip(): c for c in df.columns}

        gene_cols = [
            "gene", "genes", "gene_name", "virulence_gene", "template",
            "name", "locus", "feature", "marker", "virulence factor",
            "virulence_factor", "resistance gene", "virulence_gene_name"
        ]
        product_cols = [
            "product", "description", "function", "annotation",
            "phenotype", "locus", "feature", "database", "virulence factor",
            "virulence_factor", "predicted phenotype", "notes"
        ]
        identity_cols = [
            "identity", "%identity", "percent_identity", "identity_percent",
            "template identity", "template_identity", "identity (%)"
        ]
        coverage_cols = [
            "coverage", "%coverage", "percent_coverage", "coverage_percent",
            "template coverage", "template_coverage", "coverage (%)"
        ]
        contig_cols = [
            "sequence", "contig", "seqid", "seq_id", "chromosome", "reference",
            "reference_id", "subject", "sseqid"
        ]
        start_cols = ["start", "begin", "hit_start", "query_start", "sstart", "ref_start"]
        end_cols = ["end", "stop", "hit_end", "query_end", "send", "ref_end"]
        strand_cols = ["strand", "orientation", "direction"]

        gene_col = next((lower_cols[c] for c in gene_cols if c in lower_cols), None)
        product_col = next((lower_cols[c] for c in product_cols if c in lower_cols), None)
        identity_col = next((lower_cols[c] for c in identity_cols if c in lower_cols), None)
        coverage_col = next((lower_cols[c] for c in coverage_cols if c in lower_cols), None)
        contig_col = next((lower_cols[c] for c in contig_cols if c in lower_cols), None)
        start_col = next((lower_cols[c] for c in start_cols if c in lower_cols), None)
        end_col = next((lower_cols[c] for c in end_cols if c in lower_cols), None)
        strand_col = next((lower_cols[c] for c in strand_cols if c in lower_cols), None)

        if gene_col is None and product_col is None:
            continue

        for _, row in df.iterrows():
            gene = row.get(gene_col, "") if gene_col else ""
            product = row.get(product_col, "") if product_col else ""

            if str(gene).strip() == "" and str(product).strip() == "":
                continue

            extra_bits = []
            for c in df.columns:
                val = str(row.get(c, "")).strip()
                if val and val.lower() not in {"nan", "none"}:
                    extra_bits.append(f"{c}={val}")
            extra = "; ".join(extra_bits[:40])

            rows.append(standard_hit_row(
                sample_id=sample_id,
                species=species,
                tool=tool_name,
                database=database_name,
                gene_name=gene,
                product=product,
                contig=row.get(contig_col, "") if contig_col else "",
                start=row.get(start_col, "") if start_col else "",
                end=row.get(end_col, "") if end_col else "",
                strand=row.get(strand_col, "") if strand_col else "",
                identity=row.get(identity_col, "") if identity_col else "",
                coverage=row.get(coverage_col, "") if coverage_col else "",
                evidence=f"{tool_name} table output",
                source_file=p,
                extra=extra,
            ))

    return rows


def recursive_json_hits(obj, path=""):
    """Recursively extract candidate virulence hits from a JSON-like object."""
    hits = []

    if isinstance(obj, dict):
        keys_lower = {str(k).lower(): k for k in obj.keys()}

        gene_key = None
        for candidate in [
            "gene", "genes", "gene_name", "virulence_gene", "template",
            "name", "locus", "feature", "marker", "virulence_factor"
        ]:
            if candidate in keys_lower:
                gene_key = keys_lower[candidate]
                break

        product_key = None
        for candidate in ["product", "description", "function", "phenotype", "annotation"]:
            if candidate in keys_lower:
                product_key = keys_lower[candidate]
                break

        identity_key = None
        for candidate in ["identity", "percent_identity", "identity_percent", "template_identity"]:
            if candidate in keys_lower:
                identity_key = keys_lower[candidate]
                break

        coverage_key = None
        for candidate in ["coverage", "percent_coverage", "coverage_percent", "template_coverage"]:
            if candidate in keys_lower:
                coverage_key = keys_lower[candidate]
                break

        if gene_key is not None or product_key is not None:
            gene = obj.get(gene_key, "") if gene_key is not None else ""
            product = obj.get(product_key, "") if product_key is not None else ""

            if str(gene).strip() or str(product).strip():
                hits.append({
                    "gene": gene,
                    "product": product,
                    "identity": obj.get(identity_key, "") if identity_key is not None else "",
                    "coverage": obj.get(coverage_key, "") if coverage_key is not None else "",
                    "extra": json.dumps(obj, ensure_ascii=False)[:2000],
                    "json_path": path,
                })

        for k, v in obj.items():
            hits.extend(recursive_json_hits(v, f"{path}/{k}"))

    elif isinstance(obj, list):
        for i, item in enumerate(obj):
            hits.extend(recursive_json_hits(item, f"{path}[{i}]"))

    return hits


def parse_json_outputs(root_dir, sample_id, species, tool_name, database_name=""):
    root_dir = Path(root_dir)
    rows = []

    if not root_dir.exists():
        return rows

    files = [root_dir] if root_dir.is_file() and root_dir.suffix.lower() == ".json" else list(root_dir.rglob("*.json"))

    for p in files:
        if p.stat().st_size == 0 or p.stat().st_size > 100_000_000:
            continue

        try:
            obj = json.loads(p.read_text(encoding="utf-8", errors="replace"))
        except Exception:
            continue

        for hit in recursive_json_hits(obj):
            rows.append(standard_hit_row(
                sample_id=sample_id,
                species=species,
                tool=tool_name,
                database=database_name,
                gene_name=hit.get("gene", ""),
                product=hit.get("product", ""),
                identity=hit.get("identity", ""),
                coverage=hit.get("coverage", ""),
                evidence=f"{tool_name} JSON output",
                source_file=p,
                extra=f"json_path={hit.get('json_path', '')}; {hit.get('extra', '')}",
            ))

    return rows


# =============================================================================
# TOOL RUNNERS
# =============================================================================

def run_abricate_for_sample(sample_id, species, assembly, sample_outdir, log_dir):
    rows = []
    logs = []

    if not RUN_ABRICATE:
        return rows, logs

    if not tool_exists("abricate"):
        logs.append({
            "time": now_str(),
            "command": "abricate",
            "returncode": "",
            "stdout_file": "",
            "log_file": "",
            "status": "missing_tool",
            "message": "ABRicate executable not found in PATH.",
        })
        return rows, logs

    databases = list(ABRICATE_BASE_DATABASES)
    if is_ecoli_or_shigella_like(species):
        databases.extend(ABRICATE_ECOLI_SPECIFIC_DATABASES)

    # Avoid duplicate DBs while preserving order.
    databases = list(dict.fromkeys(databases))

    available_dbs = get_available_abricate_databases()
    if available_dbs:
        filtered = []
        for db in databases:
            if db in available_dbs:
                filtered.append(db)
            else:
                logs.append({
                    "time": now_str(),
                    "command": f"abricate --db {db}",
                    "returncode": "",
                    "stdout_file": "",
                    "log_file": "",
                    "status": "skipped_missing_database",
                    "message": f"ABRicate database '{db}' not found. Available: {', '.join(sorted(available_dbs))}",
                })
        databases = filtered

    out_files_for_summary = []

    for db in databases:
        out_file = Path(sample_outdir) / f"abricate_{db}.tsv"
        log_file = Path(log_dir) / f"{sample_id}_abricate_{db}.log"
        out_files_for_summary.append(out_file)

        cmd = [
            "abricate",
            "--db", db,
            "--minid", str(ABRICATE_MIN_IDENTITY),
            "--mincov", str(ABRICATE_MIN_COVERAGE),
            str(assembly),
        ]

        log = run_command(cmd, stdout_path=out_file, log_path=log_file)
        logs.append(log)

        if log["status"] in {"ok", "dry_run"} and not DRY_RUN_ONLY:
            rows.extend(parse_abricate_output(out_file, sample_id, species, db))

    if CREATE_ABRICATE_SUMMARY and out_files_for_summary and tool_exists("abricate") and not DRY_RUN_ONLY:
        existing_files = [str(p) for p in out_files_for_summary if Path(p).exists()]
        if existing_files:
            summary_out = Path(sample_outdir) / "abricate_summary.tsv"
            summary_log = Path(log_dir) / f"{sample_id}_abricate_summary.log"
            cmd = ["abricate", "--summary"] + existing_files
            logs.append(run_command(cmd, stdout_path=summary_out, log_path=summary_log))

    return rows, logs


def run_virulencefinder_for_sample(sample_id, species, assembly, sample_outdir, log_dir):
    rows = []
    logs = []

    if not RUN_VIRULENCEFINDER:
        return rows, logs

    # Use only for E. coli/Shigella-like genomes in this Enterobacteriaceae-focused analysis.
    if not is_ecoli_or_shigella_like(species):
        return rows, logs

    outdir = Path(sample_outdir) / "virulencefinder"
    ensure_dir(outdir)

    json_out = outdir / "virulencefinder_results.json"
    log_file = Path(log_dir) / f"{sample_id}_virulencefinder.log"

    cmd = [
        sys.executable, "-m", "virulencefinder",
        "-ifa", str(assembly),
        "-o", str(outdir),
        "-j", str(json_out),
        "-d", "all",
        "-l", str(VIRULENCEFINDER_MIN_COVERAGE),
        "-t", str(VIRULENCEFINDER_MIN_IDENTITY),
        "-x",
    ]

    if VIRULENCEFINDER_DB_PATH:
        cmd.extend(["-p", str(VIRULENCEFINDER_DB_PATH)])

    env = os.environ.copy()
    if VIRULENCEFINDER_DB_PATH:
        env["CGE_VIRULENCEFINDER_DB"] = str(VIRULENCEFINDER_DB_PATH)

    log = run_command(cmd, log_path=log_file, env=env)
    logs.append(log)

    if log["status"] in {"ok", "dry_run"} and not DRY_RUN_ONLY:
        rows.extend(parse_json_outputs(outdir, sample_id, species, "VirulenceFinder", "VirulenceFinder DB"))
        rows.extend(parse_generic_table_outputs(outdir, sample_id, species, "VirulenceFinder", "VirulenceFinder DB"))

    return rows, logs


def run_spifinder_for_sample(sample_id, species, assembly, sample_outdir, log_dir):
    rows = []
    logs = []

    if not RUN_SPIFINDER or not SPIFINDER_CMD_TEMPLATE:
        return rows, logs

    if not is_salmonella_like(species):
        return rows, logs

    outdir = Path(sample_outdir) / "spifinder"
    ensure_dir(outdir)

    command = SPIFINDER_CMD_TEMPLATE.format(
        assembly=shlex.quote(str(assembly)),
        outdir=shlex.quote(str(outdir)),
        sample=shlex.quote(str(sample_id)),
    )

    log_file = Path(log_dir) / f"{sample_id}_spifinder.log"

    log = run_command(command, log_path=log_file, shell=True)
    logs.append(log)

    if log["status"] in {"ok", "dry_run"} and not DRY_RUN_ONLY:
        rows.extend(parse_generic_table_outputs(outdir, sample_id, species, "SPIFinder", "SPIFinder DB"))

    return rows, logs


# =============================================================================
# SUMMARIES
# =============================================================================

def deduplicate_hits(df):
    if df.empty:
        return df

    tmp = df.copy()

    gene_product_key = (
        tmp["gene_name"].fillna("").astype(str).str.lower().str.replace(r"[^a-z0-9]+", "", regex=True)
        + "__"
        + tmp["product_or_locus"].fillna("").astype(str).str.lower().str.replace(r"[^a-z0-9]+", "", regex=True).str[:60]
    )

    # Once locus tags are available, deduplicate primarily by current locus_tag
    # within each sample and functional group. This is much safer for downstream
    # clustering than deduplicating only by database gene names, which can vary
    # across VFDB/VICTORS/Ecoli_VF.
    if "locus_tag" in tmp.columns:
        locus_key = tmp["locus_tag"].fillna("").astype(str).str.strip()
        tmp["dedup_key"] = locus_key.where(locus_key != "", gene_product_key)
        tmp["has_locus_tag"] = (locus_key != "").astype(int)
    else:
        tmp["dedup_key"] = gene_product_key
        tmp["has_locus_tag"] = 0

    priority = {
        "ABRicate": 1,          # ABRicate has coordinates, hence reliable locus mapping.
        "VirulenceFinder": 2,
        "SPIFinder": 3,
    }
    tmp["tool_priority"] = tmp["tool"].map(priority).fillna(9)

    tmp = tmp.sort_values(
        by=["sample_id", "functional_group", "dedup_key", "has_locus_tag",
            "tool_priority", "identity_percent", "coverage_percent"],
        ascending=[True, True, True, False, True, False, False],
    )

    dedup = tmp.drop_duplicates(
        subset=["sample_id", "species", "functional_group", "dedup_key"],
        keep="first",
    ).drop(columns=["dedup_key", "has_locus_tag", "tool_priority"], errors="ignore")

    return dedup

def make_functional_group_counts(df):
    if df.empty:
        return pd.DataFrame()

    counts = (
        df.groupby(["sample_id", "species", "functional_group"], dropna=False)["gene_name"]
        .nunique()
        .reset_index(name="unique_gene_or_locus_count")
    )

    pivot = counts.pivot_table(
        index=["sample_id", "species"],
        columns="functional_group",
        values="unique_gene_or_locus_count",
        fill_value=0,
        aggfunc="sum",
    ).reset_index()

    pivot.columns.name = None
    return pivot


def make_functional_group_gene_lists(df):
    if df.empty:
        return pd.DataFrame()

    def join_unique(values):
        vals = sorted({str(v) for v in values if str(v).strip() and str(v).lower() != "nan"})
        return "; ".join(vals)

    agg_dict = {
        "genes_or_loci": ("gene_name", join_unique),
        "products_or_loci": ("product_or_locus", join_unique),
        "tools": ("tool", join_unique),
        "databases": ("database", join_unique),
    }
    if "locus_tag" in df.columns:
        agg_dict["locus_tags"] = ("locus_tag", join_unique)
    if "EntrezGeneID" in df.columns:
        agg_dict["EntrezGeneIDs"] = ("EntrezGeneID", join_unique)

    out = (
        df.groupby(["sample_id", "species", "functional_group"], dropna=False)
        .agg(**agg_dict)
        .reset_index()
    )
    return out


def sanitize_excel_sheet_name(name: str, used_names=None) -> str:
    """Return an Excel-compatible sheet name, unique within used_names."""
    used_names = used_names if used_names is not None else set()
    safe = re.sub(r"[\\/*?:\[\]]", "_", str(name or "sheet")).strip()
    safe = safe[:31] if safe else "sheet"
    base = safe
    counter = 1
    while safe in used_names:
        suffix = f"_{counter}"
        safe = base[:31 - len(suffix)] + suffix
        counter += 1
    used_names.add(safe)
    return safe


def _unique_ordered_nonempty(values):
    seen = set()
    out = []
    for value in values:
        value = str(value or "").strip()
        if not value or value.lower() in {"nan", "none"}:
            continue
        if value not in seen:
            seen.add(value)
            out.append(value)
    return out


def make_functional_counts_workbook(dedup_hits: pd.DataFrame, output_dir: Path):
    """Create Functional counts.xlsx with one sheet per sample/pathogen.

    Each sheet contains columns of locus tags: a first 'virulence' column with all
    deduplicated virulence hits, followed by one column per functional group.
    Hits without mapped locus tags are intentionally not inserted into these
    columns, because this workbook is intended as input for locus-tag-based
    downstream clustering analyses. They remain visible in the main workbook.
    """
    if not CREATE_FUNCTIONAL_COUNTS_WORKBOOK:
        return None

    out_file = Path(output_dir) / FUNCTIONAL_COUNTS_WORKBOOK_NAME

    if dedup_hits.empty:
        with pd.ExcelWriter(out_file, engine="openpyxl") as writer:
            pd.DataFrame({"note": ["No deduplicated virulence hits were available."]}).to_excel(
                writer, sheet_name="README", index=False
            )
        return out_file

    df = dedup_hits.copy()
    if "locus_tag" not in df.columns:
        df["locus_tag"] = ""

    # Use only mapped locus tags for the clustering-ready workbook.
    df_tags = df[df["locus_tag"].fillna("").astype(str).str.strip() != ""].copy()

    used_sheet_names = set()
    try:
        writer_context = pd.ExcelWriter(out_file, engine="openpyxl")
    except PermissionError:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_file = Path(output_dir) / f"Functional counts_{timestamp}.xlsx"
        writer_context = pd.ExcelWriter(out_file, engine="openpyxl")

    with writer_context as writer:
        readme = pd.DataFrame([
            {
                "item": "Purpose",
                "description": "One sheet per sample/pathogen. Columns contain deduplicated locus tags for all virulence hits and for each functional group.",
            },
            {
                "item": "virulence column",
                "description": "All deduplicated virulence-hit locus tags for the corresponding sample.",
            },
            {
                "item": "Missing locus tags",
                "description": "Hits without mapped locus tags are excluded from these list columns and remain available in Virulence_functional_screening_results.xlsx.",
            },
            {
                "item": "Correctness",
                "description": "Locus tags are mapped only from exact annotation files using coordinate overlap. Provide GFF/GBFF files matching the exact FASTA assemblies.",
            },
        ])
        readme.to_excel(writer, sheet_name="README", index=False)
        used_sheet_names.add("README")

        missing = df[df["locus_tag"].fillna("").astype(str).str.strip() == ""].copy()
        if not missing.empty:
            keep_cols = [c for c in [
                "sample_id", "species", "gene_name", "product_or_locus", "functional_group",
                "tool", "database", "locus_mapping_status", "contig", "start", "end", "source_file"
            ] if c in missing.columns]
            missing[keep_cols].to_excel(writer, sheet_name="missing_locus_tags", index=False)
            used_sheet_names.add("missing_locus_tags")

        for sample_id, sample_df in df_tags.groupby("sample_id", sort=False):
            sample_df = sample_df.copy()
            species_values = _unique_ordered_nonempty(sample_df.get("species", []))
            species_label = species_values[0] if species_values else ""
            sheet_name = sanitize_excel_sheet_name(sample_id, used_sheet_names)

            functional_groups = _unique_ordered_nonempty(sample_df["functional_group"].tolist())
            functional_groups = sorted(functional_groups)

            columns = {
                "virulence": _unique_ordered_nonempty(sample_df["locus_tag"].tolist())
            }
            for group in functional_groups:
                group_df = sample_df[sample_df["functional_group"] == group]
                columns[group] = _unique_ordered_nonempty(group_df["locus_tag"].tolist())

            max_len = max([len(v) for v in columns.values()] + [1])
            padded = {k: v + [""] * (max_len - len(v)) for k, v in columns.items()}
            out_df = pd.DataFrame(padded)

            # Insert sample/species as first informational columns.
            out_df.insert(0, "species", [species_label] + [""] * (len(out_df) - 1))
            out_df.insert(0, "sample_id", [sample_id] + [""] * (len(out_df) - 1))
            out_df.to_excel(writer, sheet_name=sheet_name, index=False)

    return out_file


def expected_major_groups_table():
    rows = [
        {
            "species": "ecoli_o157_h7",
            "display_name": "EHEC O157:H7",
            "expected_major_groups": (
                "LEE-encoded T3SS apparatus/translocon; LEE and non-LEE T3SS effectors; "
                "intimin/Tir adhesion; Shiga toxin prophages; pO157-associated factors; "
                "fimbriae/autotransporters; LPS/O157 antigen"
            ),
            "notes": (
                "Strong T3SS + toxin profile. Use VFDB, VICTORS, Ecoli_VF, "
                "optional VirulenceFinder. Manually inspect LEE and stx prophages."
            ),
        },
        {
            "species": "ecoli_upec_uti89",
            "display_name": "UPEC UTI89",
            "expected_major_groups": (
                "Type 1 fimbriae; P fimbriae and other adhesins; siderophores/iron acquisition; "
                "toxins such as hemolysin or CNF-like factors depending on strain; serum resistance; "
                "LPS/O-antigen; capsule/surface factors"
            ),
            "notes": (
                "Use VFDB, VICTORS, Ecoli_VF, UPEC/ExPEC_VF and optional VirulenceFinder. "
                "UPEC is expected to be adhesion/iron-acquisition dominated."
            ),
        },
        {
            "species": "ecoli_upec",
            "display_name": "UPEC",
            "expected_major_groups": (
                "Fimbrial adhesins; siderophores/iron acquisition; toxins; serum resistance; "
                "LPS/O-antigen; capsule/surface factors"
            ),
            "notes": "Generic UPEC label. Prefer strain-specific labels such as ecoli_uti89 when possible.",
        },
        {
            "species": "ecoli_generic",
            "display_name": "E. coli generic",
            "expected_major_groups": (
                "Depends strongly on pathotype. Use Ecoli_VF, UPEC/ExPEC_VF and optional VirulenceFinder where possible."
            ),
            "notes": "Generic E. coli label. Manual pathotype interpretation is strongly recommended.",
        },
        {
            "species": "salmonella_typhimurium",
            "display_name": "Salmonella Typhimurium",
            "expected_major_groups": (
                "SPI-1 T3SS; SPI-2 T3SS; SPI-1/SPI-2 effectors; fimbriae; LPS; "
                "intracellular survival and stress genes; macrophage survival genes; virulence regulators"
            ),
            "notes": "Use VFDB/VICTORS plus optional SPIFinder or manual SPI curation.",
        },
        {
            "species": "salmonella_enterica",
            "display_name": "Salmonella enterica",
            "expected_major_groups": (
                "Pathogenicity-island-encoded T3SS systems depending on serovar; fimbriae; LPS; "
                "intracellular survival genes; virulence regulators"
            ),
            "notes": "Generic Salmonella label. Serovar-specific interpretation is recommended.",
        },
        {
            "species": "shigella_flexneri",
            "display_name": "Shigella flexneri",
            "expected_major_groups": (
                "pINV Mxi-Spa T3SS; Ipa translocators/effectors; Ipg chaperones; Osp/IpaH effectors; "
                "IcsA/VirG actin-based motility; LPS/O-antigen"
            ),
            "notes": "Virulence plasmid loss or fragmentation in assemblies can cause false negatives.",
        },
        {
            "species": "klebsiella_pneumoniae",
            "display_name": "Klebsiella pneumoniae",
            "expected_major_groups": (
                "Capsule/K locus; LPS/O locus; siderophores including enterobactin, yersiniabactin, "
                "aerobactin and salmochelin; rmpA/rmpA2 hypermucoidy; fim/mrk fimbriae; "
                "colibactin in some lineages; T6SS in some strains"
            ),
            "notes": "Use this ABRicate-based output as a first-pass screen. For Klebsiella-specific capsule, O-locus and acquired virulence typing, analyze the genome separately with the companion Klebsiella-specific script.",
        },
        {
            "species": "enterobacter_cloacae",
            "display_name": "Enterobacter cloacae complex",
            "expected_major_groups": (
                "Less standardized. Likely surface polysaccharides/LPS, fimbriae/adhesion, iron acquisition, "
                "T6SS and stress/host fitness genes depending on strain."
            ),
            "notes": "Species-specific virulence databases are limited; database-driven calls are incomplete and need manual/literature curation.",
        },
    ]
    return pd.DataFrame(rows)


def notes_table():
    rows = [
        {
            "topic": "Exhaustiveness",
            "note": "No database-based screen is truly exhaustive. The output is a high-confidence, database-driven candidate virulence list.",
        },
        {
            "topic": "Strain specificity",
            "note": "Virulence repertoires differ strongly by strain, plasmid content, prophages and assembly quality.",
        },
        {
            "topic": "Fragmented assemblies",
            "note": "Draft assemblies may fragment genes/loci, reducing coverage and causing false negatives.",
        },
        {
            "topic": "Functional group classification",
            "note": "Functional groups are rule-based and should be manually curated before publication.",
        },
        {
            "topic": "E. coli pathotypes",
            "note": "EHEC O157:H7 and UPEC UTI89 are both E. coli, but their virulence strategies differ strongly. Interpret pathotype-specific outputs separately.",
        },
        {
            "topic": "Enterobacter cloacae",
            "note": "Enterobacter-specific virulence databases are limited; absence of evidence is not evidence of absence.",
        },
        {
            "topic": "SPIFinder",
            "note": "SPIFinder local command-line installation varies; this script includes a command-template hook but leaves it disabled by default.",
        },
        {
            "topic": "Unsupported species",
            "note": "This version focuses on Enterobacterales/Enterobacteriaceae-like pathogens and intentionally skips Enterococcus.",
        },
    ]
    return pd.DataFrame(rows)


# =============================================================================
# METADATA LOADING
# =============================================================================

def load_metadata(metadata_file):
    metadata_file = Path(metadata_file)
    if not metadata_file.exists():
        raise FileNotFoundError(f"Metadata file not found: {metadata_file}")

    # Detect Excel files even if accidentally named .csv.
    with open(metadata_file, "rb") as f:
        magic = f.read(4)

    if metadata_file.suffix.lower() in [".xlsx", ".xlsm", ".xls"] or magic.startswith(b"PK"):
        df = pd.read_excel(metadata_file, dtype=str).fillna("")
    else:
        # sep=None lets pandas infer comma/semicolon/tab.
        df = pd.read_csv(metadata_file, dtype=str, sep=None, engine="python").fillna("")

    # Clean column names.
    df.columns = [str(c).strip() for c in df.columns]

    required = {"sample_id", "species", "assembly"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            "Metadata file is missing required columns: "
            f"{missing}. Found columns: {list(df.columns)}"
        )

    # Keep only required columns plus any extra columns for reporting.
    df["sample_id"] = df["sample_id"].apply(safe_sample_id)
    df["species_raw"] = df["species"].astype(str)
    df["species"] = df["species"].apply(normalize_species_name)
    df["assembly"] = df["assembly"].apply(lambda x: str(Path(windows_path_to_wsl(x))))

    # Optional exact annotation file used to map ABRicate coordinates to current locus tags.
    # If no annotation column exists, try to auto-detect an annotation file with the same basename.
    annotation_files = []
    for _, row in df.iterrows():
        annotation_files.append(find_annotation_file_for_assembly(row.get("assembly", ""), row))
    df["annotation_file"] = annotation_files

    # Remove empty rows.
    df = df[df["sample_id"].astype(str).str.strip() != ""].copy()
    df = df[df["assembly"].astype(str).str.strip() != ""].copy()

    # Skip unsupported species if requested.
    if SKIP_UNSUPPORTED_SPECIES:
        unsupported = df[~df["species"].isin(SUPPORTED_SPECIES)].copy()
        if not unsupported.empty:
            print("\n[WARNING] Skipping unsupported species rows:")
            print(unsupported[["sample_id", "species_raw", "species"]].to_string(index=False))
        df = df[df["species"].isin(SUPPORTED_SPECIES)].copy()

    # Validate files.
    missing_files = []
    for _, row in df.iterrows():
        if not Path(row["assembly"]).exists():
            missing_files.append((row["sample_id"], row["assembly"]))

    if missing_files:
        msg = "Assembly file(s) not found:\n"
        for sample_id, assembly in missing_files:
            msg += f"  - {sample_id}: {assembly}\n"
        msg += "\nCheck the assembly paths in input_metadata.csv and the exact FASTA filenames."
        raise FileNotFoundError(msg)

    print("\nLoaded metadata:")
    metadata_display_cols = ["sample_id", "species_raw", "species", "assembly", "annotation_file"]
    print(df[metadata_display_cols].to_string(index=False))

    no_ann = df[df["annotation_file"].astype(str).str.strip() == ""]
    if ENABLE_LOCUS_TAG_MAPPING and not no_ann.empty:
        print("\n[WARNING] No annotation file found for these samples. Locus tags/EntrezGeneID cannot be guaranteed for them:")
        print(no_ann[["sample_id", "assembly"]].to_string(index=False))
        print("Provide a matching .gff/.gff3/.gbff/.gbk file in input_metadata.csv using an 'annotation' column, or place it next to the FASTA with the same basename.")

    return df


# =============================================================================
# MAIN
# =============================================================================

def main():
    global VIRULENCEFINDER_DB_PATH

    metadata_file, output_dir, vf_db_path = resolve_runtime_paths()
    VIRULENCEFINDER_DB_PATH = vf_db_path

    output_dir = Path(output_dir)
    raw_dir = output_dir / "raw_tool_outputs"
    log_dir = output_dir / "logs"

    ensure_dir(output_dir)
    ensure_dir(raw_dir)
    ensure_dir(log_dir)

    print("\nVirulence functional screening pipeline")
    print("======================================")
    print(f"Python executable: {sys.executable}")
    print(f"Metadata file:     {metadata_file}")
    print(f"Output folder:     {output_dir}")
    if RUN_VIRULENCEFINDER and VIRULENCEFINDER_DB_PATH:
        print(f"VirulenceFinder DB: {VIRULENCEFINDER_DB_PATH}")

    metadata = load_metadata(metadata_file)

    all_rows = []
    all_logs = []

    for _, sample in metadata.iterrows():
        sample_id = sample["sample_id"]
        species = sample["species"]
        assembly = Path(sample["assembly"])

        print(f"\n[{now_str()}] Processing {sample_id} ({species})")
        sample_outdir = raw_dir / sample_id
        ensure_dir(sample_outdir)

        rows, logs = run_abricate_for_sample(sample_id, species, assembly, sample_outdir, log_dir)
        all_rows.extend(rows)
        all_logs.extend(logs)

        rows, logs = run_virulencefinder_for_sample(sample_id, species, assembly, sample_outdir, log_dir)
        all_rows.extend(rows)
        all_logs.extend(logs)

        rows, logs = run_spifinder_for_sample(sample_id, species, assembly, sample_outdir, log_dir)
        all_rows.extend(rows)
        all_logs.extend(logs)

    raw_hits = pd.DataFrame(all_rows)
    logs_df = pd.DataFrame(all_logs)

    if raw_hits.empty:
        print("\n[WARNING] No virulence hits were parsed. Check tool logs and database installation.")
        raw_hits = pd.DataFrame(columns=[
            "sample_id", "species", "tool", "database", "gene_name", "product_or_locus",
            "functional_group", "contig", "start", "end", "strand", "identity_percent",
            "coverage_percent", "accession", "evidence", "source_file", "extra",
            "locus_tag", "EntrezGeneID", "locus_mapping_status", "locus_mapping_method",
            "annotation_file", "mapped_feature_type", "mapped_gene_name", "mapped_product",
            "mapped_feature_start", "mapped_feature_end", "mapped_feature_strand",
            "overlap_bp", "overlap_fraction_of_hit", "overlap_fraction_of_feature"
        ])
    else:
        raw_hits = annotate_hits_with_locus_tags(raw_hits, metadata)

    dedup_hits = deduplicate_hits(raw_hits) if WRITE_DEDUPLICATED_TABLE else raw_hits.copy()
    counts = make_functional_group_counts(dedup_hits)
    gene_lists = make_functional_group_gene_lists(dedup_hits)
    expected = expected_major_groups_table()
    notes = notes_table()

    excel_out = output_dir / "Virulence_functional_screening_results.xlsx"

    try:
        with pd.ExcelWriter(excel_out, engine="openpyxl") as writer:
            raw_hits.to_excel(writer, sheet_name="all_raw_hits", index=False)
            dedup_hits.to_excel(writer, sheet_name="deduplicated_hits", index=False)
            counts.to_excel(writer, sheet_name="functional_counts", index=False)
            gene_lists.to_excel(writer, sheet_name="functional_gene_lists", index=False)
            metadata.to_excel(writer, sheet_name="input_metadata", index=False)
            logs_df.to_excel(writer, sheet_name="tool_run_log", index=False)
            expected.to_excel(writer, sheet_name="expected_major_groups", index=False)
            notes.to_excel(writer, sheet_name="notes", index=False)
    except PermissionError:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        excel_out = output_dir / f"Virulence_functional_screening_results_{timestamp}.xlsx"
        with pd.ExcelWriter(excel_out, engine="openpyxl") as writer:
            raw_hits.to_excel(writer, sheet_name="all_raw_hits", index=False)
            dedup_hits.to_excel(writer, sheet_name="deduplicated_hits", index=False)
            counts.to_excel(writer, sheet_name="functional_counts", index=False)
            gene_lists.to_excel(writer, sheet_name="functional_gene_lists", index=False)
            metadata.to_excel(writer, sheet_name="input_metadata", index=False)
            logs_df.to_excel(writer, sheet_name="tool_run_log", index=False)
            expected.to_excel(writer, sheet_name="expected_major_groups", index=False)
            notes.to_excel(writer, sheet_name="notes", index=False)
        print("\n[WARNING] Original Excel output was open/locked. Wrote timestamped file instead.")

    csv_out = output_dir / "Virulence_functional_screening_deduplicated_hits.tsv"
    dedup_hits.to_csv(csv_out, sep="\t", index=False)

    functional_counts_out = make_functional_counts_workbook(dedup_hits, output_dir)

    print("\nDone.")
    print(f"Raw hits:       {len(raw_hits)}")
    print(f"Dedup hits:     {len(dedup_hits)}")
    if "locus_tag" in dedup_hits.columns:
        n_locus = int((dedup_hits["locus_tag"].fillna("").astype(str).str.strip() != "").sum())
        print(f"Dedup hits with locus tags: {n_locus}")
    print(f"Excel output:   {excel_out}")
    print(f"TSV output:     {csv_out}")
    if functional_counts_out is not None:
        print(f"Functional counts workbook: {functional_counts_out}")
    print(f"Logs:           {log_dir}")
    print("\nPlease inspect the 'tool_run_log' sheet first to identify any failed optional tool calls.")
    print("For reliable locus tags, inspect 'locus_mapping_status' in all_raw_hits/deduplicated_hits.")


if __name__ == "__main__":
    main()
