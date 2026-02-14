# ProteoMapper

ProteoMapper is a Python-based desktop application for analyzing protein sequence alignments from Excel files.
It provides an intuitive graphical interface to clean sequences, expand alignments, detect motifs, highlight positions, scan Pfam domains using HMMER, and quantify motif–domain relationships.

The tool is designed for researchers working with protein families, conserved motifs, and domain annotations who want fast visual inspection and structured Excel outputs.

---

## ✨ Features

* 📄 Reads Excel files with automatic header detection
* 🧬 Cleans protein sequences (handles `>` headers, gaps, invalid characters)
* 🔢 Expands sequences into aligned per-position columns in Excel
* 🔍 Regex-based motif detection with highlighting
* 📊 Match summary with counts and positions
* 📌 Position-based residue highlighting
* 🧠 Optional conservation threshold for frequent motifs
* 🧩 Pfam domain scanning using **HMMER (`hmmscan`)**
* 📈 Domain summary and annotated alignment sheet
* 🧮 **MDCS summary** for motif–domain overlap interpretation
* ⚡ Parallel processing for faster motif and domain analysis
* 🖥️ Responsive Tkinter GUI with background execution

---

## 📦 Requirements

### Python

* Python **3.8+**

### Core packages

```
pandas
numpy
openpyxl
```

Tkinter is bundled with most Python installations.

### Optional (for domain scanning)

* **HMMER** installed
* `hmmscan` available on PATH
* Pfam database file (`.hmm`)

---

## 🚀 Installation

Clone repository:

```
git clone <your-repo-url>
cd proteomapper
```

Install dependencies:

```
pip install -r requirements.txt
```

---

## ▶️ Quick Start

Run the GUI:

```
python proteomapper/ProteoMapper.py
```

Then:

1. Select your Excel input file
2. Choose analysis options
3. Click **Proceed**
4. Output Excel will be generated next to your input file

---

## 📑 Input Excel Format

ProteoMapper scans the first sheet and automatically detects the header row.

Required columns:

* `Protein Sequences`
* Either `Gene Name` or `Gene ID`

Example:

```
Gene Name | Protein Sequences
GeneA     | MKT--AILVGL
GeneB     | MKTAKAILVGL
```

Notes:

* Leading `>` is automatically removed
* Gaps (`-` or spaces) are preserved in alignment
* Non-amino-acid characters are removed

---

## 🔎 Analysis Options

### Motif Search (Regex)

Enter one regex per line in the GUI.

Example:

```
C..C
N[^P][ST]
```

Results:

* Matches highlighted in Excel
* Counts and positions listed in **Match Summary**

You can also define a **conservation threshold (%)** to outline frequently occurring motifs.

---

### Position Highlighting

Enter space-separated residue positions (1-based):

```
5 10 25 42
```

These columns will be highlighted across all sequences.

---

### Domain Scanning (Pfam / HMMER)

ProteoMapper can run:

```
hmmscan Pfam-A.hmm sequences.fasta
```

Outputs:

* Domain-highlighted alignment sheet
* Domain summary table
* Annotated Excel comments for each domain

---

## 🧠 MDCS Summary (Motif–Domain Context Score)

ProteoMapper computes MDCS values describing how motifs overlap with domains:

* **1.0 → Fully embedded in domain**
* **0 < MDCS < 1 → Partial overlap**
* **0 → Outside domains**
* **No Domains → sequence has no detected domains**

This helps identify:

* domain-specific motifs
* conserved functional regions
* motif–domain relationships in protein families

---

## 🛠️ Installing HMMER

### Check installation

```
hmmscan -h
```

### Conda (recommended)

```
conda install -c bioconda hmmer
```

### macOS

```
brew install hmmer
```

### Ubuntu/Debian

```
sudo apt-get install hmmer
```

---

## 🗂 Outputs

Generated next to your input file:

* `msa_data_processed.xlsx` → main output

### Key worksheets

* **First_Sheet** → cleaned alignment + highlights
* **Match Summary** → motif counts and ranges
* **Domain_Highlights** → alignment with domain annotation
* **Domain Summary** → detailed domain hits
* **MDCS Summary** → motif–domain overlap interpretation

---

## 🧪 Testing

Run:

```
pytest -q
```

Some tests may require:

* GUI environment
* HMMER installed

---

## 🧯 Troubleshooting

**Header not detected**

* Ensure column names match:

  * `Protein Sequences`
  * `Gene Name` or `Gene ID`

**HMMER not found**

* Ensure `hmmscan` is installed and on PATH

**Slow performance**

* Reduce motif patterns
* Increase parallel processes


## License
This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

Copyright (c) 2025 Sifullah Mahmud Sefa.

Made with ❤️ for the bioinformatics community.


