# TAD Benchmarking

This repository contains code for benchmarking and analyzing Topologically Associating Domains (TADs) in genomic data.
It includes implementations of various metrics for evaluating TAD callers and analyzing their characteristics.

## Project Overview

TADs (Topologically Associating Domains) are self-interacting genomic regions that play important roles in gene
regulation and 3D genome organization. This project provides a comprehensive benchmarking framework for evaluating
different TAD calling algorithms and methods.

The analysis focuses on several key aspects of TAD performance:

1. **Histone Modification Enrichment (Fig1C)** - Evaluating TAD with histone modifications (H3K27me3, H3K36me3)
2. **Structural Protein Enrichment (Fig1D)** - Analyzing enrichment of structural proteins (e.g., CTCF, RAD21, SMC3) at
   TAD boundaries
3. **CTCF Orientation Analysis (Fig1E)** - Examining CTCF orientation at TAD boundaries
4. **Contact Enrichment (Fig1F)** - Measuring interaction frequencies within TADs
5. **Boundary Insulation (Fig1G)** - Quantifying the insulation strength of TAD boundaries
6. **Robustness (Fig1H)** - Assessing the consistency of TAD calls across different Hi-C sequencing depths

## Repository Structure

- `Data/` - Contains some example data
- `Fig1C_HistMod_FDR/` - Histone modification analysis
- `Fig1D_StructProt_EnrichBoundaries/` - Structural protein enrichment at TAD boundaries
- `Fig1E_CTCF_orientation/` - CTCF orientation analysis
- `Fig1F_Contact_Enrichment/` - Contact enrichment analysis
- `Fig1G_Boundary_Insulation/` - Boundary insulation analysis
- `Fig1H_Robustness/` - Robustness analysis
- `utils/` - Common utility functions used across analyses

## Requirements

The code is written in R and requires the following packages:

- magrittr
- dplyr
- tibble
- readr
- tidyr
- purrr
- data.table
- GenomicRanges
- Biostrings
- BSgenome.Hsapiens.UCSC.hg19
- TFBSTools
- plyr
- stringr

## Installation

```r
# Install CRAN packages
install.packages(c("magrittr", "dplyr", "tibble", "readr", "tidyr", "purrr",
                   "data.table", "plyr", "stringr"))

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GenomicRanges", "Biostrings",
                       "BSgenome.Hsapiens.UCSC.hg19", "TFBSTools"))
```

## Usage

Each analysis component can be run independently by executing the main R script in the corresponding directory.

### Example: Running Histone Modification Analysis

```r
source("Fig1C_HistMod_FDR/Fig1C_HistMod_FDR.R")
```

### Example: Running CTCF Orientation Analysis

```r
source("Fig1E_CTCF_orientation/Fig1E_CTCF_orientation.R")
```

### Note on Hi-C Matrix Data

For **Contact Enrichment (Fig1F)** and **Boundary Insulation (Fig1G)** analyses, Hi-C contact matrix files are required
but not included in this repository due to their large size.

Please download the Hi-C matrix data from:
**TADShop: [https://tadshop.unil.ch/](https://tadshop.unil.ch/)**

The expected file format is: `HiC_matrix_{chr}.txt` (e.g., `HiC_matrix_chr6.txt`)

## Input Data Format

### TAD Data (`Data/TADs.tsv`)

The TAD data is stored in TSV format with the following columns:

- `chr` - Chromosome name (e.g., chr1, chr2, ..., chrX)
- `start` - Start position of the TAD (in base pairs, 1-based)
- `end` - End position of the TAD (in base pairs)
- `meta.tool` - Name of the TAD caller tool used
- `meta.resol` - Resolution of the TAD calls (in base pairs)
- `meta.sub` - Subsampling identifier (used in robustness analysis; "." indicates full dataset)

## Contributing

Contributions to improve the benchmarking framework are welcome. Please feel free to submit pull requests or open issues
for discussion.

## Citation

If you use this benchmarking framework in your research, please cite our work.

## License

This project is licensed under the terms of the included LICENSE file.

## Contact

For questions or issues, please open an issue on GitHub or visit [TADShop](https://tadshop.unil.ch/).
