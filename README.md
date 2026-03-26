┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃                                                    ProteoMapper                                                    ┃
┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛

ProteoMapper is a Python-based desktop application for integrated analysis of protein sequence alignments, motif
detection, and domain annotation. It provides an intuitive graphical interface to process protein sequences, detect
regulatory motifs, scan Pfam domains using HMMER, and quantify motif–domain spatial relationships with alignment-aware
visualization.

The tool is designed for researchers working with protein families, conserved motifs, and domain annotations who need
reproducible, alignment-centric analysis with structured Excel outputs.

──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

                                                   ✨ Key Features

                                                Input & Preprocessing

 • 📄 Dual input format support: Excel (.xlsx, .xls) and FASTA (.fasta, .fa, .faa)
 • 🔍 Automatic header detection for flexible Excel formats
 • 🧬 Intelligent sequence cleaning: handles gaps, invalid characters, and diverse formats
 • 🔢 Alignment expansion: per-position columns for visual inspection

                                                    Motif Analysis

 • 🎯 Regex-based motif detection with full pattern control
 • 📊 Positional conservation scoring: identify evolutionarily constrained motifs
 • 📈 Positional dispersion (σ): quantify motif placement variability
 • 🧠 Conservation threshold: highlight motifs above user-defined frequency
 • 🔖 Match summary: counts, positions, and conservation statistics

                                                  Domain Annotation

 • 🧩 Pfam domain scanning using HMMER (hmmscan)
 • 🏷️ Interactive domain visualization with embedded metadata
 • 📋 Domain summary tables with E-values and bit scores
 • ⚡ Parallel processing for large-scale analyses

                                               Motif–Domain Integration

 • 🧮 MDCS (Motif–Domain Coverage Score): quantify spatial overlap
 • 📊 Three-table MDCS summary:
    • Per-sequence motif–domain relationships
    • Domain-specific motif enrichment patterns
    • Motif-level statistics (mean/median MDCS, SLiM identification)
 • 🎯 Distinguish domain-embedded motifs from SLiMs in disordered regions

                                                   User Experience

 • 🖥️ Responsive Tkinter GUI with progress tracking
 • 📦 Unified Excel output with multiple analysis sheets
 • ⚡ Multiprocessing support for faster computation
 • 📌 Custom position highlighting for functional sites

──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

                                                   📦 Requirements

                                                        Python

 • Python 3.8+ (3.9+ recommended)

                                                  Core Dependencies

```
pandas>=1.3.0
numpy>=1.21.0
openpyxl>=3.0.0
```

Tkinter is bundled with most Python installations

                                            Optional (for domain scanning)

 • HMMER 3.x (hmmscan accessible via PATH)
 • Pfam-A database (.hmm file with pressed indices)

──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

                                                   🚀 Installation

                                                   Clone Repository

```bash
git clone https://github.com/yourusername/proteomapper.git
cd proteomapper
```

                                                 Install Dependencies

```bash
pip install -r requirements.txt
```

                                               Verify HMMER (optional)

```bash
hmmscan -h
```

──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

                                                    ▶️ Quick Start

                                                      Launch GUI

```bash
python ProteoMapper.py
```

                                                    Basic Workflow

 1 Select input file (Excel or FASTA)
 2 Configure analysis options:
    • Add motif patterns (optional)
    • Set conservation threshold (default: 60%)
    • Enable domain scanning (optional)
 3 Click Proceed
 4 Output Excel generated in same directory as input

──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

                                                   📑 Input Formats

                                                     Excel Format

ProteoMapper automatically detects header rows containing:

Required columns:

 • Protein Sequences (or Sequences)
 • At least one identifier: Gene Name, Gene ID, Protein Name, or Protein ID

Example:

```
Gene Name | Protein Sequences
GeneA     | MKT--AILVGL
GeneB     | MKTAKAILVGL
```

Notes:

 • Gaps (-) preserved for alignment visualization
 • Non-amino-acid characters automatically removed
 • Multiple identifier columns supported

──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
                                                     FASTA Format

Standard FASTA files with flexible header parsing:

Supported formats:

```
>sp|P12345|PROT1 Protein description
>P12345 Protein description
>P12345
```

Parsing behavior:

 • First field or pipe-delimited ID → Protein ID
 • Remaining text → Protein Name
 • Handles UniProt, NCBI, and custom formats

──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

                                                 🔎 Analysis Options

                                               1. Motif Search (Regex)

Define patterns in GUI (one per line):

Examples:

```
C..C              # Zinc finger motif
N[^P][ST][^P]     # N-glycosylation site
[ST]..E           # Acidic phosphorylation motif
...([ST])P..      # Proline-directed kinase site
```

Outputs:

 • Sky-blue highlighting in alignment
 • Match Summary table with counts and positions
 • Positional dispersion (σ) for each pattern

──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
                                              2. Conservation Threshold

Set minimum percentage (default: 60%) for positional conservation.

Interpretation:

 • Motifs at ≥ threshold → Red border (positionally conserved)
 • σ = 0 → Strict positional constraint
 • σ > 0 → Positional variability despite sequence conservation

──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
                                               3. Position Highlighting

Highlight specific alignment columns (space-separated, 1-based):

```
5 10 25 42 67
```

Displays as green fill across all sequences.

──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
                                              4. Domain Scanning (HMMER)

Requirements:

 • HMMER installed (hmmscan on PATH)
 • Pfam-A.hmm database with indices

Execution:

```bash
hmmscan --domtblout output.txt Pfam-A.hmm sequences.fasta
```

Outputs:

 • Orange-highlighted domains in alignment
 • Hoverable comments with metadata (E-value, bit score, accession)
 • Domain Summary table

──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

                                    🧮 MDCS Summary (Motif–Domain Coverage Score)

Definition:

```
MDCS = (motif–domain overlap length) / (motif length)
```

Range: 0 to 1

 • 1.0 → Fully embedded in domain
 • 0 < MDCS < 1 → Partial overlap
 • 0 → Outside all domains

                                             Three-Table Output Structure

Table 1: Per-Sequence MDCS

 • Sequence ID, Motif Pattern, MDCS values, Interpretation

Table 2: Domain-Specific Motif Embedding

 • Domain, Motif Pattern, % Full/Partial embedding, Count
 • Identifies motif-rich domains

Table 3: Motif-Level Statistics

 • Motif Pattern, Mean MDCS, Median MDCS, Count
 • Outside Domains count (SLiM indicator)
 • Touching Domains count

Biological Interpretation:

 • Low mean/median MDCS → Likely SLiMs in disordered regions
 • High mean/median MDCS → Domain-associated signatures
 • Complements positional conservation analysis

──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

                                                  🗂 Output Structure

                                                    Generated File

<input_filename>_processed.xlsx (Excel workbook)

                                                      Worksheets

                                                                                                         
  Sheet               Description                                                                        
 ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 
  First_Sheet         Aligned sequences with motif highlights (sky-blue) and conservation borders (red)  
  Match_Summary       Motif counts, positions, conservation %, positional dispersion (σ)                 
  Domain_Highlights   Alignment with orange domain annotations and metadata comments                     
  Domain_Summary      E-values, bit scores, coordinates for all detected domains                         
  MDCS_Summary        Three tables: per-sequence, domain-enrichment, motif statistics                    
                                                                                                         

                                                   Additional Files

 • hitdata.txt → Raw HMMER domain table output

──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

                                                  🛠️ Installing HMMER

                                                  Check Installation

```bash
hmmscan -h
```

                                                 Conda (Recommended)

```bash
conda install -c bioconda hmmer
```

                                                  Platform-Specific

macOS (Homebrew):

```bash
brew install hmmer
```

Ubuntu/Debian:

```bash
sudo apt-get install hmmer
```

Windows: Download precompiled binaries from hmmer.org

                                                    Pfam Database

Download from: ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/

```bash
# Press database
hmmpress Pfam-A.hmm
```

──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

                                                      🧪 Testing

```bash
# Run all tests
pytest

# Quick test
pytest -q

# With coverage
pytest --cov=ProteoMapper
```

Test suite includes:

 • FASTA parsing validation
 • Excel input processing
 • Motif detection accuracy
 • Domain scanning integration
 • MDCS calculation correctness

──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

                                                  🧯 Troubleshooting

                                                "Header not detected"

Solution: Ensure column names match:

 • Protein Sequences (case-insensitive)
 • At least one: Gene Name, Gene ID, Protein Name, Protein ID

                                                  "HMMER not found"

Solution:

```bash
which hmmscan  # Should return path
export PATH=$PATH:/path/to/hmmer/bin
```

                                          "Maximum columns exceeded" (Excel)

Issue: Sequences > 16,384 amino acids exceed Excel limit

Solutions:

 • Extract region of interest
 • Use CSV output (future feature)
 • Split alignment into chunks

                                                   Slow Performance

Optimizations:

 • Reduce regex patterns
 • Increase parallel processes (GUI option)
 • Use E-value cutoff for domain scanning (default: 0.001)

                                                 FASTA parsing issues

Common fixes:

 • Ensure headers start with >
 • Check for special characters in headers
 • Verify sequences contain only amino acids and gaps

──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

                                                     📚 Citation

If you use ProteoMapper in your research, please cite:

```
[Your manuscript citation will go here]
```

──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

                                                   🤝 Contributing

Contributions are welcome! Please:

 1 Fork the repository
 2 Create a feature branch (git checkout -b feature/YourFeature)
 3 Commit changes (git commit -m 'Add YourFeature')
 4 Push to branch (git push origin feature/YourFeature)
 5 Open a Pull Request

──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────


## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

Copyright (c) 2025 Sifullah Mahmud Sefa.

Made with ❤️ for the bioinformatics community.


