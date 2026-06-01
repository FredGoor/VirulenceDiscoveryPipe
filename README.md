[README.md](https://github.com/user-attachments/files/28460361/README.md)
# Virulence Discovery Pipeline

Python pipeline for identifying candidate virulence genes and virulence-associated functional groups in assembled bacterial genomes.

The repository contains two complementary scripts:

1. **`GF11_abricate_virulence_functional_screen_publishable.py`**  
   ABRicate-based virulence gene screening across Enterobacterales/Enterobacteriaceae-like pathogens.

2. **`GF12_kleborate_virulence_gene_extractor_publishable.py`**  
   Kleborate-guided extraction of *Klebsiella pneumoniae* virulence-associated locus tags from a matching genome annotation file.

## Short description

This pipeline screens assembled bacterial genomes for virulence-associated genes and loci, classifies hits into broad biological functions, and exports publication-ready Excel/TSV tables for downstream analyses such as functional clustering, codon-usage analyses, or comparative virulence profiling.

## Main features

- Runs ABRicate on assembled bacterial genomes using virulence databases such as VFDB, VICTORS, Ecoli_VF and UPEC_ExPEC_VF.
- Supports metadata-driven batch analysis of multiple genomes.
- Optionally parses VirulenceFinder and SPIFinder outputs when these tools are installed and configured locally.
- Maps ABRicate hits to current locus tags using matching GFF/GFF3/GBK/GBFF annotation files.
- Deduplicates hits across databases and tools.
- Assigns hits to broad functional virulence categories.
- Generates Excel workbooks and TSV files suitable for manual curation and downstream analysis.
- Provides a separate Kleborate-guided workflow for *Klebsiella pneumoniae* species complex genomes.

## Recommended workflow

1. Assemble genomes and collect matching annotation files.
2. Prepare `input_metadata.csv`.
3. Run the ABRicate-based screen with GF11.
4. Inspect the `tool_run_log` sheet for failed tool calls or missing databases.
5. Inspect `locus_mapping_status` to verify locus-tag mapping.
6. Manually curate important categories and ambiguous hits.
7. For *Klebsiella pneumoniae* species complex genomes, run GF12 separately.
8. Use the deduplicated tables and functional locus-tag lists for downstream analyses.

### Python packages

The scripts require Python 3.9 or later and the following Python packages:

```bash
conda create -n virulence_screen python=3.11 pandas openpyxl biopython
conda activate virulence_screen
```

### External command-line tools

Depending on which scripts and options are used, install:

- ABRicate
- ABRicate virulence databases such as VFDB, VICTORS, Ecoli_VF and UPEC_ExPEC_VF
- Kleborate
- VirulenceFinder, optional
- SPIFinder, optional

The external tools should be available in the active environment, meaning commands such as the following should work from the terminal:

```bash
abricate --version
abricate --list
kleborate --version
```

## Input files

### GF11 ABRicate-based screen

The main input is a metadata table with at least three columns:

```csv
sample_id,species,assembly
Ecoli_O157H7,ecoli_o157_h7,genomes/Ecoli_O157H7.fna
SL1344,salmonella_typhimurium,genomes/SL1344.fna
Sflex,shigella_flexneri,genomes/Shigella_flexneri.fna
```

To enable locus-tag mapping, add a matching annotation file:

```csv
sample_id,species,assembly,annotation
Ecoli_O157H7,ecoli_o157_h7,genomes/Ecoli_O157H7.fna,annotations/Ecoli_O157H7.gff
SL1344,salmonella_typhimurium,genomes/SL1344.fna,annotations/SL1344.gff
Sflex,shigella_flexneri,genomes/Shigella_flexneri.fna,annotations/Shigella_flexneri.gff
```

Accepted metadata formats:

- `.csv`
- `.tsv`
- semicolon-separated `.csv`
- `.xlsx`

Accepted annotation formats:

- `.gff`
- `.gff3`
- `.gbk`
- `.gbff`
- `.gb`
- `.genbank`

### GF12 Kleborate-guided Klebsiella extractor

Required inputs:

1. *Klebsiella pneumoniae* species complex assembly FASTA file.
2. Matching annotation file in GFF/GFF3/GBK/GBFF format.
3. Output folder.

## Usage

### Interactive mode with file pickers

Both scripts can be launched without arguments. File and folder pickers will open.

```bash
conda activate virulence_screen
python GF11_abricate_virulence_functional_screen_publishable.py
```

```bash
conda activate virulence_screen
python GF12_kleborate_virulence_gene_extractor_publishable.py
```

## GF11 outputs

The ABRicate-based script writes the following files to the selected output folder:

```text
Virulence_functional_screening_results.xlsx
Virulence_functional_screening_deduplicated_hits.tsv
Functional_counts.xlsx
raw_tool_outputs/
logs/
```

The main Excel workbook contains:

| Sheet | Description |
|---|---|
| `all_raw_hits` | All parsed hits from all tools/databases before deduplication. |
| `deduplicated_hits` | Deduplicated virulence hits, preferably resolved to locus tags when annotation files are available. |
| `functional_counts` | Count matrix of functional virulence groups per sample. |
| `functional_gene_lists` | Gene/locus lists grouped by sample and functional category. |
| `input_metadata` | Cleaned metadata table used by the script. |
| `tool_run_log` | Command lines, return codes and error messages for each external tool call. |
| `expected_major_groups` | Reference table of expected major virulence groups by pathogen/pathotype. |
| `notes` | Caveats and interpretation notes. |

The additional `Functional_counts.xlsx` workbook contains one sheet per sample, with columns of locus tags suitable for downstream clustering or enrichment analyses.

## GF12 outputs

The Kleborate-guided script writes:

```text
Kleborate_extracted_Klebsiella_virulence_genes.xlsx
kleborate_run.log
```

The Excel workbook contains:

| Sheet | Description |
|---|---|
| `Kleborate_raw` | Raw parsed Kleborate table. |
| `Kleborate_modules` | Detected/reported Kleborate virulence modules. |
| `Extracted_genes` | Extracted genes and locus tags from the annotation file. |
| `Summary_counts` | Number of extracted locus tags per category and confidence class. |
| `Locus_tag_lists` | Category-wise locus-tag lists. |
| `Notes` | Interpretation notes and relevant caveats. |
| category-specific sheets | One sheet per virulence category. |

## Functional categories

The scripts use rule-based functional grouping. Categories include:

- Type III secretion apparatus/translocon
- Type III secretion effectors
- Type VI secretion / interbacterial competition
- Iron acquisition / siderophores
- Capsule / surface polysaccharide
- LPS / O-antigen
- Adhesins / fimbriae / pili
- Biofilm-associated factors
- Motility / cell-to-cell spread
- Intracellular survival / stress resistance
- Immune evasion / serum resistance
- Virulence regulation
- Toxins and genotoxins
- Pathogenicity islands / virulence plasmids
- Unknown virulence-associated hits

These categories are designed for first-pass interpretation and should be manually curated before final biological conclusions are drawn.


