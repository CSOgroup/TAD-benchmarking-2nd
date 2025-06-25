# TAD Benchmarking

This repository contains code for benchmarking and analyzing Topologically Associating Domains (TADs) in genomic data.
It includes implementations of various metrics for evaluating TAD callers and analyzing their characteristics.

## Project Overview

TADs (Topologically Associating Domains) are self-interacting genomic regions that play important roles in gene
regulation and 3D genome organization. This project provides a comprehensive benchmarking framework for evaluating
different TAD calling algorithms and methods.

The analysis focuses on several key aspects of TAD performance:

1. **Histone Modification Enrichment** - Evaluating TAD with histone modifications (H3K27me3, H3K36me3)
2. **Structural Protein Enrichment** - Analyzing enrichment of structural proteins (e.g., CTCF, RAD21, SMC3) at TAD
   boundaries
3. **Contact Enrichment** - Measuring interaction frequencies within TADs
4. **Boundary Insulation** - Quantifying the insulation strength of TAD boundaries
5. **Robustness** - Assessing the consistency of TAD calls across different Hi-C sequencing depth

## Repository Structure

- `Data/` - Contains some example data
- `Fig1C_HistMod_FDR/` - Histone modification analysis
- `Fig1D_StructProt_EnrichBoundaries/` - Structural protein enrichment at TAD boundaries
- `Fig1E_Contact_Enrichment/` - Contact enrichment analysis
- `Fig1F_Boundary_Insulation/` - Boundary insulation analysis
- `Fig1G_Robustness/` - Robustness analysis
- `utils/` - Common utility functions used across analyses

## Requirements

The code is written in R and requires the following packages:

- magrittr
- dplyr
- tibble
- readr
- GenomicRanges
- purrr

## Usage

Each analysis component can be run independently by executing the main R script in the corresponding directory.

## Input Data Format

The TAD data is stored in TSV format with the following columns:

- `chr` - Chromosome name
- `start` - Start position of the TAD
- `end` - End position of the TAD
- `meta.tool` - Name of the TAD caller tool used
- `meta.resol` - Resolution of the TAD calls (in base pairs)
- `meta.sub` - Used in robustness analysis

## Contributing

Contributions to improve the benchmarking framework are welcome. Please feel free to submit pull requests or open issues
for discussion.

## License

This project is licensed under the terms of the included LICENSE file. 