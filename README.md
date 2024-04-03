## Overview

This repository contains codes for collecting global genome assemblies of Escherichia coli from diverse sources, annotating these assemblies with the blaCTX-M-15 gene, and subsequently classifying the contigs as either chromosomal or plasmid.

## Workflow

1. **Data Scraping**: Genome assemblies of *E. coli* are scraped from diverse sources, encompassing public databases, research articles, and sequencing projects, to form a global collection of strains.

2. **Metadata Curation**: Relevant metadata associated with each genome assembly, including geographic location, host species, isolation source, and antibiotic resistance profiles, are curated and integrated into the analysis, providing essential context for understanding *E. coli* genetic diversity and epidemiology.

3. **Annotation**: The blaCTX-M-15 gene is annotated within each genome assembly using Abricate (Resfinder database) for sequence alignment and homology detection. This annotation offers insights into the presence or absence of this antibiotic resistance gene across different strains.

4. **Contig Classification**: Following annotation, contigs within each assembly are classified as either chromosomal or plasmid based on sequence features, such as size, GC content, and presence of plasmid-associated genes. This classification aids in understanding the genomic context and potential mobility of the blaCTX-M-15 gene.

## Code Structure

- `modules/`: Contains the processes for the Nexflow subworkflows.
- `conda/`: Contains yaml files for the conda environments.
- `bin/`: Contains individual scripts for parsing and plotting all relevant data.

