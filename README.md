# Gene-Level SNP Grouping, Insertion Mapping & Heritability Toolkit

![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue)

A small toolkit for:

- **Extracting** SNPs for arbitrary gene regions from large HapMap/VCF files  
- **Collapsing** diploid calls to haploid representations and dominant genotypes  
- **Grouping** strains by identical SNP‐match profiles (network-based clustering)  
- **Visualising** representative strain sets and match tables as heat-map PNGs  and circos plotting


The code was developed for *Bacillus subtilis* genome‐wide association analyses but is file-format-agnostic and should work with any species as long as your HapMap/VCF files follow standard conventions.

---

---

## Repository layout
```
.
├── conftest.py
├── heritability_illustrator
│   ├── __init__.py
│   ├── class_analysis.py             # Class/phenotype group analysis helpers
│   ├── convert_vcf.py                # Core VCF parsing and normalization
│   ├── create_gene_region_profile.py # SNP extraction for gene/region intervals
│   ├── create_visualization.py       # Heatmap + profile visualization pipeline
│   ├── data_loaders.py               # File I/O utilities for VCF and metadata
│   ├── df_helpers.py                 # Pandas dataframe cleaning / merging
│   ├── dominant_matching.py          # Collapse genotypes to dominant matches
│   ├── haploid_squashing.py          # Convert diploid calls to haploid reps
│   ├── heterozygote_to_homozygote.py # Handle heterozygous calls simplification
│   ├── load_and_extract_data.py      # High-level data extraction orchestration
│   ├── network_helpers.py            # Graph/network SNP clustering functions
│   ├── plot_images.py                # Matplotlib/seaborn visualization wrappers
│   └── utils.py                      # Small utility functions
├── notebooks
│   └── example_notebook.ipynb        # Walkthrough notebook
├── README.md
└── requirements.txt```


---

## Installation

### 1. Clone and create an environment

```bash

git clone https://github.com/LaneMatthewJ/heritability_illustrator.git
cd heritability_illustrator
python -m venv .venv        # or conda create -n genekit python=3.11
source .venv/bin/activate
```

or with conda: 

```bash
conda create -n genekit python=3.11
conda activate genekit
```

### 2. Istall Python Requirements: 
```bash
pip install -r requirements.txt
```

**Note**: There may be a potential optional dependency: 
```bash
sudo apt-get install wkhtmltopdf
```

## Example Usage

```python
import heritability_illustrator
from heritability_illustrator.create_visualization import generate_profile_from_vcf, create_start_stop_df

start_pos = 0
stop_pos = 4_200_000

outfiles_ronn1_ME = generate_profile_from_vcf(
    "../variant_mapping/Bsubtilis_350strain_SNPs_final.vcf", 
    start_pos = start_pos,
    stop_pos =  stop_pos,
    class_encoding_filename = None, # '../image_data/sorted/class_encodings.xlsx'
    input_image_directory = None) 

create_start_stop_df(outfiles_ronn1_ME[0], start_pos, stop_pos, "Bsubtilis_SNPs", 0) 

df = pd.read_csv("Bsubtilis_350strain_SNPs_final_0_to_4200000_visualizations/Bsubtilis_350strain_SNPs_final_0_4200000_snp_matchings_by_group.tsv", sep='\t', index_col=0) 
data = df[df.columns[:-1]].loc[df.index[1:]]
colors = ['cyan', 'black']

black_to_blue_cmap = LinearSegmentedColormap.from_list('black_to_blue', colors)
sns.heatmap(data, cmap=black_to_blue_cmap) 
```

##  Example Applications
-	Mapping gene-level SNP differences between experimental strains.
-	Identifying strain clusters with identical genomic variation patterns.
-	Performing heritability-informed grouping of phenotypes vs. genotypes.
-	Generating publication-ready heatmaps and circos plots.
