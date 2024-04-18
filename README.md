## Overview

This repository contains codes for collecting global genome assemblies of Escherichia coli from diverse sources, annotating these assemblies with the blaCTX-M-15 gene, and subsequently classifying the contigs as either chromosomal or plasmid and investigating their phylogenetic relationships.

## Code Structure

- `modules/`: Contains the processes for the Nexflow subworkflows.
- `conda/`: Contains yaml files for the conda environments.
- `bin/`: Contains individual scripts for parsing and plotting all relevant data.

## Workflow

1. **Data Scraping**: Genome assemblies of *E. coli* are scraped from diverse sources, encompassing public databases, research articles, and sequencing projects, to form a global collection of strains.

2. **Metadata Curation**: Relevant metadata associated with each genome assembly, including geographic location, host species, isolation source, and antibiotic resistance profiles, are curated and integrated into the analysis, providing essential context for understanding *E. coli* genetic diversity and epidemiology.

  Fetch metadata form NCBI for the Blackwell samples
```
nextflow run -entry fetch_ncbimetadata --inputCsv N5037_ctm15pos_Blackwell.tsv --outdir ./
```
or
```
python3 bin/run_esummary.py
```
4. **Annotation**: The blaCTX-M-15 gene is annotated within each genome assembly using Abricate (Resfinder database) for sequence alignment and homology detection. This annotation offers insights into the presence or absence of this antibiotic resistance gene across different strains.

5. **Contig Classification**: Following annotation, contigs within each assembly are classified as either chromosomal or plasmid based on sequence features, such as size, GC content, and presence of plasmid-associated genes. This classification aids in understanding the genomic context and potential mobility of the blaCTX-M-15 gene.

6. **Flanking Sequence Extraction**: Flanking sequences spanning 5 kilobases upstream and downstream of the blaCTX-M-15 gene are extracted from the genome assemblies. These sequences provide context for understanding the genetic environment surrounding the gene and potential regulatory elements or mobile genetic elements associated with its dissemination.

7. **Dating Phylogenies**: A phylogenetic tree was constructed using the extracted flanking sequences from *E. coli* genomes. Bayesian phylogenetic inference methods are employed to estimate the divergence times between different lineages, providing insights into the evolutionary history and temporal dynamics of *E. coli* populations. This analysis aids in understanding the spread and emergence of antibiotic resistance within the context of *E. coli* evolution.

  
  Extract individual gene products from Genbank files
```
python3 bin/extract_gene_product.py $file.gbk -p "IS1380 family transposase ISEcp1" "Beta-lactamase CTX-M-1"
```
Batch run:
```
bash bin/batch_extract_gene_product.sh
```
  Extract segments from ISEcp1 to blaCTX-M-15 and collapse under one header
```
python3 bin/extract_startgene_endgene_collapse.py $file.gbk -s "IS1380 family transposase ISEcp1" -e "Beta-lactamase CTX-M-1"
```
Don't collapse:
```
python3 bin/extract_startgene_endgene.py $file.gbk -s "IS1380 family transposase ISEcp1" -e "Beta-lactamase CTX-M-1"
```
Batch run:
```
bash bin/batch_extract_StartEnd_genes.sh
```
or
```
nextflow run --entry extractgenes --genbank_dir /path/to/genbank/files --outdir /path/to/output/directory 
```
Concatenate the extracted segments for the next steps e.g. progressive alignment and phylogeny generation
```
cat *_cut.fasta > concatenated_flank_cut.fasta
```

