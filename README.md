# Gene-Level SNP Grouping, Insertion Mapping & Heritability Toolkit

![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue)
![License: MIT](https://img.shields.io/badge/license-MIT-green)

A small but opinionated toolkit for:

* **Extracting** SNPs for arbitrary gene regions from large HapMap/VCF files  
* **Collapsing** diploid calls to haploid representations and dominant genotypes  
* **Grouping** strains by identical SNP‐match profiles (network-based clustering)  
* **Visualising** representative strain sets and match tables as heat-map PNGs  and circos plotting


The code was developed for *Bacillus subtilis* genome‐wide association analyses but is file-format-agnostic and should work with any species as long as your HapMap/VCF files follow standard conventions.

---

## Table of Contents
1. [Repository layout](#repository-layout)  
2. [Installation](#installation)  
3. [Quick start](#quick-start)  
4. [Input file formats](#input-file-formats)  
5. [Outputs](#outputs)  
6. [Library API](#library-api)  
7. [Notebook workflow](#notebook-workflow-heritability_estimationipynb)  
8. [Troubleshooting](#troubleshooting)  
9. [Contributing](#contributing)  
10. [License & citation](#license--citation)

---

## Repository layout
```
├── MAIN.py                         # Entry-point helper script
├── heritability_estimation.ipynb   # End-to-end heritability notebook
├── lib/                            # Core library (importable)
│   ├── init.py
│   ├── load_and_extract_data.py
│   ├── dominant_matching.py
│   ├── create_gene_region_profile.py
│   ├── plot_images.py
│   └── … (helpers: data_loaders, df_helpers, etc.)
├── data/                           # Example HapMap / VCF files 
└── README.md
```

> **Tip** The library does **not** depend on the notebook; all heavy-lifting lives in `lib/`.

---

## Installation

### 1. Clone and create an environment

```bash
git clone https://github.com/<your-org>/gene-snp-toolkit.git
cd gene-snp-toolkit
python -m venv .venv        # or conda create -n genekit python=3.11
source .venv/bin/activate
```
### 2. Istall Python Requirements: 
```bash
pip install -r requirements.txt
```

### 3. External Dependencies: 
```bash
sudo apt-get install wkhtmltopdf
```





## Quick Start: 

Below is a minimal end-to-end example run directly from Python. Replace file paths and genomic coordinates with your own.

```python
start_pos = 0
stop_pos = 4_200_000
proportion_similarity_for_clustering=0

profile = generate_profile_from_vcf(
    "./data/Bsubtilis_350strain_SNPs_final.vcf", 
    start_pos = start_pos,
    stop_pos =  stop_pos
)

create_start_stop_df(
  profile[0], 
  start_pos, 
  stop_pos, 
  "Bsubtilis_SNPs", 
  proportion_similarity_for_clustering
) 
```

All generated artefacts are collected in a folder named

`<gene>_<start>_to_<stop>_visualizations/`





