# Automated assembly of molecular mechanisms at scale from text mining and curated databases

This repository contains scripts, notebooks, and other files
for generating the results of the manuscript.

The INDRA software itself is available at: https://github.com/sorgerlab/indra.

## Reading, databases, and assembly

We created a Benchmark Corpus of Statements using INDRA. As input for text
mining, a list of publications was obtained as follows:

* reading/pmids.txt
  * Obtained in 2019 by querying Entrez for PubMed publications on each of a
    list of human genes to get a set of 567507 PMIDs.
  * Shuffle the PMIDs in a list and save the entries into pmids.txt.

The output of reading systems as well as the content of multiple curated
databases is processed in the `run_assembly` module. This module also runs
a set of assembly steps that filter and transform the raw Statements into
non-redundant assembled Statements.

The `run_assembly` folder has its own Makefile specifically for performing
the assembly.

A top-level `Makefile` automates the workflow for generating various statistics and figure panels based on the assembled Benchmark Corpus Statements.

## Curation and belief models

The `bioexp/curation` module contains functions to manage the curated
corpus and train various belief models on the curation data.

## Analysis of the Benchmark Corpus

The `bioexp/figures` module contains scripts to generate figure panels for
Figures 2, 3, and 4.

Multiple notebooks perform analysis on the Benchmark Corpus to produce
metrics and figures used in the manuscript:
* `notebooks/Statement distribution.ipynb`: analysis on the distribution of mentions per Statement for Figure 2.
* `notebooks/Reader_overlap_and_error_analysis.ipynb`: analysis based on curation of multi-reader overlap and empirical correctness for Figure 5.
* `notebooks/Biogrid benchmark.ipynb`: analysis for the "Validation of assembled mechanisms and comparison against curated resources" section and Figure 6.
* `notebooks/Depmap Benchmark.ipynb`: analysis for the "Detecting and explaining gene dependency correlations with an assembled causal network" section and Figure 7. The notebook also makes use of code in the `bioexp/depmap` module to preprocess the DepMap data.

## Data files

* The Benchmark Corpus is available as a Python pickle of INDRA Statement objects at https://bioexp-paper.s3.amazonaws.com/bioexp_asmb_preassembled.pkl.
* Aggregated curations on the Benchmark Corpus are available as a Python pickle file at https://bioexp-paper.s3.amazonaws.com/curation_dataset_with_bg_psp.pkl. The content is a list, where each element of the list is a dictionary prociding a `stmt_hash` corresponding to the hash of a Statement in the Benchmark Corpus, and giving a `correct` flag (0 for overall incorrect, 1 for overall correct) to represent the overall result of curation, along with some other metadata.
