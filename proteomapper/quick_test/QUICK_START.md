# ProteoMapper — Quick Start / Installation Check

This folder contains a minimal, synthetic dataset for verifying that a local
ProteoMapper installation is working correctly, before running the tool on a
real dataset. It takes under a minute to run.

## Contents

| File | Purpose |
|---|---|
| `sample_alignment.fasta` | 6 short aligned toy protein sequences (FASTA input) |
| `sample_alignment.xlsx` | The same 6 sequences (Excel input) |
| `example_patterns.txt` | 3 example regex motif patterns, one per line |
| `expected_output.xlsx` | The verified correct output for this input, for comparison |
| `QUICK_START.md` | This file |

The dataset is intentionally tiny and has no biological significance — it
exists only to exercise the motif-matching and conservation-highlighting
logic on a fast, human-checkable example. It does **not** cover domain
scanning (HMMER/Pfam), which depends on a local Pfam database and is not
practical to bundle here; see the main repository README for domain-scanning
setup instructions.

## Steps

1. Launch ProteoMapper.
2. Under **1. Select Input File**, click **Browse** and select either
   `sample_alignment.fasta` or `sample_alignment.xlsx` (both give identical
   results — this is a good way to confirm both input parsers are working).
3. Under **2. Select Analysis Options**, check **Motif Searching** only.
   Leave **Position Highlighting** and **Domain Scanning** unchecked.
4. In the **Motif Search Options** panel that appears:
   - Open `example_patterns.txt` in a text editor, copy its contents, and
     paste them into the **Regex Patterns** box (one pattern per line).
   - Leave **Conservation Threshold (%)** at its default value of `60`.
5. Click **Proceed**. Processing should complete in a few seconds.
6. Open the generated `msa_data_processed.xlsx` and compare it against
   `expected_output.xlsx` in this folder.

## What a correct run looks like

- **`First_Sheet`**: a 6-row alignment matrix (columns 1–20), with:
  - Columns **7–9** (`KRS`) highlighted sky-blue in 5 of 6 sequences, with a
    thick red border — this pattern (`K[RK]S`) is deliberately absent from
    `GENE5`, and is conserved at ≥60% of sequences.
  - Columns **12–13** (`PG`) highlighted sky-blue with a thick red border in
    **all** 6 sequences.
  - No highlighting on `GENE5` at columns 7–9, confirming the mismatch is
    handled correctly.
- **`Match Summary`** sheet listing all three patterns:

  | Regex Pattern | Matches | Positional Dispersion (σ) |
  |---|---|---|
  | `K[RK]S` | 5 | 0.00 |
  | `PG` | 6 | 0.00 |
  | `AA` | 3 | 1.41 |

  The `AA` pattern is included specifically to demonstrate a motif that
  falls **below** the 60% conservation threshold (so it is highlighted
  sky-blue but receives no red border) and that occurs at more than one
  alignment position (hence a nonzero positional dispersion, σ = 1.41).

If your output matches this pattern of highlighting and these Match Summary
values, your installation is working correctly and you can proceed to your
own dataset.

## If something doesn't match

- **No output file is created / an error dialog appears:** check that all
  required Python packages are installed and that you are running the
  correct entry point for the application.
- **Motif columns are not highlighted:** confirm the patterns were pasted
  into the Regex Patterns box exactly as they appear in
  `example_patterns.txt` (one pattern per line, including the `::Label`
  suffix), and that **Motif Searching** is checked.
- **Different Match counts than shown above:** confirm you loaded
  `sample_alignment.fasta` or `sample_alignment.xlsx` from this folder
  unmodified, and that the Conservation Threshold field was left at `60`.

For domain-scanning (HMMER/Pfam) verification, or for questions beyond this
quick check, see the main repository README and the manuscript
(Section 2.3–2.4) for full setup and usage details.
