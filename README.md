# Butyrate_Propionate_Comparative_Analysis
This repository contains the reproducible bioinformatics pipelines used in the manuscript:

“Analysis of butyrate and propionate pathway abundance in human gut microbiota reveals family-level separation”
doi: https://doi.org/10.1101/2023.11.27.568948

The workflows implemented here identify butyrate- and propionate-producing pathways in bacterial genomes and estimate their relative abundance in metagenomic samples using profile HMMs and read mapping–based normalization.

The full methodological rationale and validation are described in the manuscript and Supplementary Information; this repository provides the exact computational implementation used for the analyses.

Overview of the Pipeline

As outlined in Supplementary Fig. S1 of the manuscript, the analysis consists of five major steps:

1. Identification of pathway-positive genomes (IMG)
Bacterial genomes from the Integrated Microbial Genomes (IMG) database are screened for the presence of required KEGG Orthology (KO) genes defining eight butyrate and propionate fermentation pathways.

2. Construction and application of HMMs (HMMER)
Profile Hidden Markov Models (HMMs) are built for each pathway gene and used to identify pathway-positive strains from a curated set of gut and reference genomes.

3. Gene catalog construction
Pathway-positive genomes are used to generate a non-redundant gene catalog containing >2,000 pathway gene sequences.

4. Read mapping to metagenomes (Bowtie2)
Metagenomic reads are mapped to the gene catalog to quantify pathway gene hits across samples.

5. Pathway abundance estimation and normalization
Pathway gene counts are normalized by gene length, pathway length, and the housekeeping gene rplB to estimate the fraction of genomes encoding each pathway in a metagenomic sample.

Steps 2–5 are implemented as Snakemake workflows with per-pipeline dependencies managed through isolated conda environments; Step 1 is a standalone R-based preprocessing step.

Execution Order

The pipelines are designed to be run in the following order:
1. IMG genome search
2. HMMER pathway identification
3. Bowtie2 read mapping and abundance estimation

Each pipeline is self-contained, with its own Snakemake workflow and dependency specifications. Although each pipeline can be executed independently, downstream pipelines assume the outputs of upstream steps.

Requirements
Software:
  - Snakemake ≥ 7
  - Conda / Mamba (recommended for environment resolution)
  - Linux or macOS
Snakemake will automatically create all required environments from the YAML files located in each pipeline’s workflow/envs/ directory.

Installation

Clone the repository:
```
git clone https://github.com/RebbyAnne/Butyrate_Propionate_Comparative_Analysis.git
cd Butyrate_Propionate_Comparative_Analysis
```
  
Install Snakemake (example using conda):
```
conda create -n snakemake -c conda-forge snakemake
conda activate snakemake
```

Running the Pipelines
1. IMG Genome Search Pipeline

This directory implements Step 1 of the analysis workflow: identification and curation of bacterial genomes from the Integrated Microbial Genomes (IMG) database that encode complete butyrate or propionate fermentation pathways.Unlike the downstream analysis steps, this component is not a Snakemake workflow. It is an R-based data filtering and curation pipeline designed to operate on genome and gene tables exported directly from the IMG database.

Purpose

The goal of this step is to define a high-confidence set of pathway-positive genomes and gene sequences that are subsequently used to:
  - Construct profile Hidden Markov Models (pHMMs)
  - Generate the gene catalog used for metagenomic read mapping
Pathway presence is determined using biologically informed completeness criteria derived from experimentally characterized model strains and manual curation of gene annotations, as described in the manuscript Methods.

Inputs
Input files are expected to reside in the img_genome_with_gene_ids/ directory and include, for each pathway gene:
  - IMG genome identifiers
  - IMG gene identifiers
  - Gene descriptions

Processing
The R Markdown script img_genome_filter.Rmd performs the following operations:
  - Loads IMG genome–gene tables for butyrate and propionate pathway genes
  - Applies pathway-specific gene completeness requirements
  - Implements manual annotation filtering for ambiguous gene descriptions where necessary
  - Excludes genes or subunits that were empirically shown to cause misclassification
  - Identifies genomes encoding complete fermentation pathways
  - Extracts curated IMG gene ID lists for each pathway gene

Distinct filtering logic is applied for each propionate (WWC, SP, Pdu) and butyrate pathway in accordance with the criteria described in the manuscript.

Outputs

Filtered gene ID lists are written to the img_gene_ids/ directory. Each output file contains IMG gene identifiers corresponding to a specific pathway gene and serves as direct input for downstream pHMM construction and HMMER searches.

Execution

This step is intended to be run once prior to executing the Snakemake-based pipelines.

The script can be executed by opening img_genome_filter.Rmd in RStudio and running all code chunks sequentially, or by rendering it via the command line using rmarkdown::render(). No Snakemake configuration is required for this stage.

2. HMMER Search Pipeline

This step constructs profile HMMs for pathway genes and identifies butyrate- and propionate-producing strains.
```
cd ../hmmer_search_pipeline
snakemake --use-conda --cores 8
```

Key outputs:
  - HMM profiles for each pathway gene
  - Filtered lists of pathway-positive strains
  - Gene catalogs derived from validated pathway-positive genomes
  - HMMER hit filtering follows the criteria described in the manuscript:
  - Score cutoff set to 50% of the lowest-scoring model strain
  - Resolution of multi-pathway strains based on relative HMM scores

3. Bowtie2 Read Mapping Pipeline
This step maps metagenomic reads to the pathway gene catalog and computes pathway abundance estimates.
```
cd ../bowtie_search_pipeline
snakemake --use-conda --cores 8
```

Bowtie2 is run using the --very-sensitive preset in end-to-end mode.
Outputs include:
  - Per-gene hit counts
  - Pathway-level abundance estimates
  - Normalized pathway prevalence values

Input Data
Input files (e.g., genome FASTA files, metagenomic reads) are specified within each pipeline’s Snakemake configuration and are not included in this repository due to size and licensing constraints.
See the manuscript and Supplementary Tables for accession numbers and data sources.


