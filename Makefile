BUILD := build
DATA := data
NET := networks
FIG1 := bioexp/figures/figure1
FIG2 := bioexp/figures/figure2

all: preprocessing fig2

preprocessing: \
        $(BUILD)/prior_genes.txt \
        $(BUILD)/pc_multidigraph.pkl

clean:
	cd $(BUILD); rm -rf *

fig1: $(BUILD)/fig1_pc_egfr_mapk1_paths.txt

fig2: $(BUILD)/fig2_evidence_distribution.pdf

# DATA -----------------------------------------------------------------------

$(DATA)/%:
	python -m bioexp.transfer_s3 get "$@" $(DATA)

$(DATA)/PathwayCommons9.All.hgnc.txt:
	wget -P $(DATA) http://www.pathwaycommons.org/archives/PC2/v9/PathwayCommons9.All.hgnc.txt.gz
	gunzip $@

# PREPROCESSING --------------------------------------------------------------

# The list of prior genes from the data and related sources
$(BUILD)/prior_genes.txt: process_data.py \
                          $(DATA)/Korkut\ et\ al.\ Data\ 05122017.xlsx \
                          $(DATA)/ras_pathway_proteins.csv \
                          $(DATA)/drug_grounding.csv
	python process_data.py

# Pathway Commons network parsed from the extended SIF and pickled as
# an instance of a networkx MultiDiGraph
$(BUILD)/pc_multidigraph.pkl: $(DATA)/PathwayCommons9.All.hgnc.txt \
                              $(FIG1)/find_paths.py \
                              $(BUILD)/prior_genes.txt
	python $(FIG1)/find_paths.py parse_pc

# FIGURE 1 -------------------------------------------------------------------

# Example pathfinding output over Pathway Commons
$(BUILD)/fig1_pc_egfr_mapk1_paths.txt: $(BUILD)/pc_multidigraph.pkl \
                                       $(FIG1)/find_paths.py
	python $(FIG1)/find_paths.py find_paths

# FIGURE 2 -------------------------------------------------------------------

$(BUILD)/fig2_evidence_distribution.pdf: \
        $(DATA)/bioexp_preassembled.pkl \
        $(FIG2)/preassembly_stats.py
	python -m bioexp.figures.figure2.preassembly_stats

