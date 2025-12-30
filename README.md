# ProteoMapper

ProteoMapper is a Python/Tkinter desktop app for analyzing protein sequence alignments from Excel files. It cleans and expands sequences into per-position columns, highlights motifs and positions, summarizes matches, and can optionally run Pfam domain scanning via HMMER.

## Features
- Reads Excel files with flexible header placement.
- Cleans and normalizes protein sequences (handles leading ">" and gap characters).
- Splits sequences into aligned, per-position columns for inspection in Excel.
- Motif search with regex highlighting and a match summary sheet.
- Optional conservation threshold to outline frequently matched motifs.
- Position-based highlighting for user-specified residue indices.
- Optional Pfam domain scanning via `hmmscan` with domain highlights and a hit report.
- Parallel processing for motif search and domain scanning to improve performance.

## Requirements
- Python 3.8+
- Core Python packages:
  - pandas
  - numpy
  - openpyxl
  - psutil
- Tkinter (typically bundled with Python on Windows/macOS)
- Optional: HMMER (for domain scanning)
  - `hmmscan` must be on PATH
  - Pfam database file (HMM format)

Install dependencies:

```bash
pip install -r requirements.txt
```

For development:

```bash
pip install -r requirements-dev.txt
```

## Quick Start
Run the GUI:

```bash
python proteomapper/ProteoMapper.py
```

From the GUI, select your input Excel file and choose the analysis options you need.

## Input Excel Format
ProteoMapper scans the first sheet and detects the header row automatically. Your file must include:
- A column named `Protein Sequences`
- Either `Gene Name` or `Gene ID`

Notes:
- The tool tolerates a leading `>` in sequences and will strip it.
- Gaps (`-` or spaces) are preserved in the output alignment view.
- Non-amino-acid characters are removed during cleaning.

Example header row:

```
Gene Name | Protein Sequences
```

Example rows:

```
GeneA | MKT--AILVGL
GeneB | MKTAKAILVGL
```

## Analysis Options
### Motif Search (Regex)
Provide one regex per line in the GUI. Matches are highlighted in the output and summarized.

Example patterns:

```
C..C
N[^P][ST]
```

A conservation threshold can be set to outline motif regions found in a high percentage of sequences.

### Position Highlighting
Enter space-separated positions (1-based) to highlight columns in the output alignment.

Example:

```
5 10 25 42
```

### Domain Scanning (Pfam / HMMER)
Enable domain scanning to run `hmmscan` against a Pfam database. Results are filtered by
E-value threshold and written to `hitdata.txt`. Detected domains are highlighted in Excel.

Requirements:
- `hmmscan` available on PATH
- Pfam database file path

## HMMER Installation Guide
Verify installation:

```bash
hmmscan -h
```

### Recommended: Conda (cross-platform)

```bash
conda install -c bioconda hmmer
```

### macOS (Homebrew)

```bash
brew install hmmer
```

### Ubuntu/Debian

```bash
sudo apt-get update
sudo apt-get install hmmer
```

### Pfam Database Setup
Download a Pfam HMM database (for example `Pfam-A.hmm`) and prepare it:

```bash
hmmpress Pfam-A.hmm
```

Then point the GUI to the `Pfam-A.hmm` file.

## Outputs
The tool writes outputs next to the input file:
- `msa_data_processed.xlsx`: main output with alignment, highlights, and summary sheet
- `hitdata.txt`: domain scan results (only when domain scanning is enabled)

Key worksheets:
- `First_Sheet`: cleaned and expanded alignment with highlights
- `Match Summary`: per-pattern counts and matched position ranges

## Examples
### Motif Highlight Example
Input sequences:

```
GeneX | ACDEFGHIKLMNPQRSTVWY
GeneY | ACDEYGHIKLMNPQRSTVWY
```

Regex pattern:

```
DEFG
```

Result:
- The matching region is highlighted in both sequences.
- The summary sheet lists the pattern and count of matches.

### Domain Scan Example
If `hmmscan` identifies a domain match for sequence 3 at positions 45-120, the output will:
- Highlight those positions in the alignment sheet.
- Add a line to `hitdata.txt` with sequence ID, start/end, accession, and E-value.

## Project Structure
- `proteomapper/ProteoMapper.py`: main application and analysis logic
- `proteomapper/tests/`: test suite
- `requirements.txt`: runtime dependencies
- `requirements-dev.txt`: development/test dependencies

## Testing
Run tests from the repository root:

```bash
pytest -q
```

Some tests may require a GUI-capable environment or external tools (HMMER) depending on
your platform and configuration.

## Troubleshooting
- "Could not find the correct header row": verify `Protein Sequences` and `Gene Name`/`Gene ID` headers.
- HMMER errors: confirm `hmmscan` is installed and available on PATH.
- Large datasets: increase available CPU cores or reduce motif patterns to improve performance.

## Notes
- The GUI remains responsive during processing; analysis runs in a background thread.
- Output formatting is optimized for Excel readability and downstream inspection.

## License
This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

Copyright (c) 2025 Sifullah Mahmud Sefa
