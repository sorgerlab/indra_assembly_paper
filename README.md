# bioexp_paper

## Data Sources

* Data
  * ras_pathway_proteins.csv. Derived from Frank McCormick's Ras Pathway v2.0,
    available at
    https://www.cancer.gov/research/key-initiatives/ras/ras-central/blog/2015/ras-pathway-v2
  * drug_grounding.csv. Uniprot IDs for each drug, from Korkut et al.
  * Korkut et al. Data 05122017.xlsx. Obtained from paper (URL?)

* Networks
  * PathwayCommons.All.hgnc.txt. Obtained from
    http://www.pathwaycommons.org/archives/PC2/v9/PathwayCommons.All.hgnc.txt.gz,
    on 1/8/2018.

## Reading

* reading/pmids.txt
  * Run indra/models/hgnc_all/get_pmids.py (commit 771f50) to get
     pmids_for_gene.pkl based on 41150 genes 
  * Collect all PMIDs into a set to get 567507 PMIDs
  * Shuffle the PMIDs in a list and save the entries into pmids.txt

## Derived Data Files

* prior_genes.txt. Produced by running process_data.py

