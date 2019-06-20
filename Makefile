BUILD := build
DATA := data
NET := networks
FIG1 := bioexp/figures/figure1
FIG2 := bioexp/figures/figure2
FIG4 := bioexp/figures/figure4
DEPLOY := ../bioexp_manuscript/figures/figure_panels

all: preprocessing fig2 fig4 fig5

sample: $(BUILD)/reach_sample_uncurated.tsv

preprocessing: \
        $(BUILD)/prior_genes.txt \
        $(BUILD)/pc_multidigraph.pkl

deploy:
	rsync -av $(BUILD)/*.pdf $(DEPLOY)

clean:
	cd $(BUILD); rm -rf *

fig1: $(BUILD)/fig1_pc_egfr_mapk1_paths.txt

fig2: $(BUILD)/fig2_evidence_distribution.pdf \
      $(BUILD)/fig2_stmt_counts_before_pa.pdf

fig4: $(BUILD)/fig4_belief_surface.pdf

fig5: $(BUILD)/reach_complexes_raw.csv

# MAKEFILE GRAPH
graph: makegraph.pdf

makegraph.pdf: makegraph.dot
	dot -T pdf makegraph.dot -o makegraph.pdf

makegraph.dot: Makefile
	make -Bnd | make2graph > makegraph.dot


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


# STMT SAMPLE FOR CURATION --------------------------------------------------
$(BUILD)/reach_sample_uncurated.tsv: $(DATA)/bioexp_filter_top_level.pkl
	python -m bioexp.curation.sample $< 100 $(BUILD) \
                             reach sparser rlimsp isi medscan trips


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

$(BUILD)/fig2_stmt_counts_before_pa.pdf: \
        $(DATA)/bioexp_preassembled.pkl \
        $(DATA)/fig2_stmt_counts.txt \
        $(FIG2)/plot_stmt_counts.py
	python -m bioexp.figures.figure2.plot_stmt_counts

# FIGURE 4 -------------------------------------------------------------------
$(BUILD)/fig4_belief_surface.pdf: \
        $(FIG4)/belief_surface.py
	python -m bioexp.figures.figure4.belief_surface

# FIGURE 5 -------------------------------------------------------------------

$(BUILD)/reach_complexes_raw.csv: $(DATA)/bioexp_reach.pkl
	python -m bioexp.figures.figure5.validate_complex \
                             $(DATA)/bioexp_reach.pkl \
                             $(BUILD)/reach_complexes_raw.tsv


