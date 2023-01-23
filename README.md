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

* Fig 4B and Table 2: run `make output/fig4_reach_curve.pdf` (which runs
  `bioexp/figures/figure4/curated_correctness.py`)
* Figure 4D and Table 3: run `make output/fig4_{reader}_model_fits.pdf` for each
  reader.
  * The makefile first runs `bioexp/curation/process_curations.py` for the
    reader to generate fits, then dumps the maximum likelihood results table
    (Table 3).
  * Then the Makefile runs `bioexp/figures/figure4/model_fit_plots.py` to
    generate the plots in Fig 4D.
  * Figures EV3A, EV3B, and EV3C (binomial, betabinomial and belief evidence
    fits): `bioexp/curation/process_curations.py`

Multiple notebooks perform analysis on the Benchmark Corpus to produce
metrics and figures used in the manuscript:
* `notebooks/Statement distribution.ipynb`: analysis on the distribution of mentions per Statement for Figure 2.
* `notebooks/Reader_overlap_and_error_analysis.ipynb`: analysis based on curation of multi-reader overlap and empirical correctness for Figure 5.
* `notebooks
* `notebooks/Biogrid benchmark.ipynb`: analysis for the "Validation of assembled mechanisms and comparison against curated resources" section and Figure 6.
* `notebooks/Depmap Benchmark.ipynb`: analysis for the "Detecting and explaining gene dependency correlations with an assembled causal network" section and Figure 7. The notebook also makes use of code in the `bioexp/depmap` module to preprocess the DepMap data.


## Data files

The Benchmark Corpus is available on Zenodo at https://zenodo.org/record/7559353.
* The Benchmark Corpus on Zenodo is called `indra_benchmark_corpus.pkl`. For reproducing results, put `indra_benchmark_corpus.pkl` in the `data`
folder and rename it to `bioexp_asmb_preassembled.pkl`.
* Curations on the Benchmark Corpus are available as JSON file on Zenodo
and in this repository: [indra_assembly_curations.json](https://github.com/sorgerlab/indra_assembly_paper/blob/master/data/curation/indra_assembly_curations.json). These raw curations are processed and aggregated to create two pickle files that are used in
the various notebooks. These files are in version control in this repository as well: [multireader_curation_dataset.pkl](https://github.com/sorgerlab/indra_assembly_paper/blob/master/data/curation/multireader_curation_dataset.pkl) and [extended_curation_dataset.pkl](https://github.com/sorgerlab/indra_assembly_paper/blob/master/data/curation/extended_curation_dataset.pkl). 
