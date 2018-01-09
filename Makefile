BUILD := build
DATA := data
NET := networks
FIG1 := figures/figure1

clean:
	cd $(BUILD); rm -rf *

all: fig1

fig1: $(BUILD)/fig1_pc_egfr_mapk1_paths.txt


# The list of prior genes from the data and related sources
$(BUILD)/prior_genes.txt: process_data.py \
                          $(DATA)/Korkut\ et\ al.\ Data\ 05122017.xlsx \
                          $(DATA)/ras_pathway_proteins.csv \
                          $(DATA)/drug_grounding.csv
	python process_data.py

# Pathway Commons network parsed from the extended SIF and pickled as
# an instance of a networkx MultiDiGraph
$(BUILD)/pc_multidigraph.pkl: $(NET)/PathwayCommons9.All.hgnc.txt \
                              $(FIG1)/find_paths.py \
                              $(BUILD)/prior_genes.txt
	python $(FIG1)/find_paths.py parse_pc


# Example pathfinding output over Pathway Commons
$(BUILD)/fig1_pc_egfr_mapk1_paths.txt: $(BUILD)/pc_multidigraph.pkl \
                                       $(FIG1)/find_paths.py
	python $(FIG1)/find_paths.py find_paths

