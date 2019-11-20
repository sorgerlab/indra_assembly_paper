# Must match config.json
OUTPUT := output
DATA := data
NET := networks
FIG1 := bioexp/figures/figure1
FIG2 := bioexp/figures/figure2
FIG4 := bioexp/figures/figure4
DEPLOY := ../bioexp_manuscript/figures/figure_panels

all: preprocessing fig2 fig4 fig5 korkut_pysb

sample: $(OUTPUT)/reach_sample_uncurated.tsv

preprocessing: \
        $(OUTPUT)/prior_genes.txt \
        $(OUTPUT)/pc_multidigraph.pkl

deploy:
	rsync -av $(OUTPUT)/*.pdf $(DEPLOY)

clean:
	cd $(OUTPUT); rm -rf *

fig1: $(OUTPUT)/fig1_pc_egfr_mapk1_paths.txt

fig2: $(OUTPUT)/fig2_evidence_distribution.pdf \
      $(OUTPUT)/fig2_stmt_counts_before_pa.pdf

fig4: $(OUTPUT)/fig4_belief_surface.pdf

fig5: $(OUTPUT)/reach_complexes_raw.csv

korkut_pysb: $(OUTPUT)/bioexp_test_paths.pkl

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
#
# The list of prior genes from the data and related sources
$(OUTPUT)/prior_genes.txt: $(DATA)/Korkut\ et\ al.\ Data\ 05122017.xlsx \
                          $(DATA)/ras_pathway_proteins.csv \
                          $(DATA)/drug_grounding.csv
	python -m bioexp.explanation.process_data \
        $(DATA)/Korkut\ et\ al.\ Data\ 05122017.xlsx \
        $(DATA)/ras_pathway_proteins.csv \
        $(DATA)/drug_grounding.csv \
        $@

# Pathway Commons network parsed from the extended SIF and pickled as # an instance of a networkx MultiDiGraph
$(OUTPUT)/pc_multidigraph.pkl: $(DATA)/PathwayCommons9.All.hgnc.txt \
                              $(FIG1)/find_paths.py \
                              $(OUTPUT)/prior_genes.txt
	python $(FIG1)/find_paths.py parse_pc


# STMT SAMPLE FOR CURATION --------------------------------------------------
$(OUTPUT)/reach_sample_uncurated.tsv: $(DATA)/bioexp_asmb_preassembled.pkl
	python -m bioexp.curation.sample $< 40 1 10 $(OUTPUT) \
                             reach sparser rlimsp isi medscan trips

# FIGURE 1 -------------------------------------------------------------------

# Example pathfinding output over Pathway Commons
$(OUTPUT)/fig1_pc_egfr_mapk1_paths.txt: $(OUTPUT)/pc_multidigraph.pkl \
                                       $(FIG1)/find_paths.py
	python $(FIG1)/find_paths.py find_paths

# FIGURE 2 -------------------------------------------------------------------

$(OUTPUT)/fig2_evidence_distribution.pdf: \
        $(DATA)/bioexp_preassembled.pkl \
        $(FIG2)/preassembly_stats.py
	python -m bioexp.figures.figure2.preassembly_stats

$(OUTPUT)/fig2_stmt_counts_before_pa.pdf: \
        $(DATA)/bioexp_preassembled.pkl \
        $(DATA)/fig2_stmt_counts.txt \
        $(FIG2)/plot_stmt_counts.py
	python -m bioexp.figures.figure2.plot_stmt_counts

# FIGURE 4 -------------------------------------------------------------------
$(OUTPUT)/fig4_belief_surface.pdf: \
        $(FIG4)/belief_surface.py
	python -m bioexp.figures.figure4.belief_surface

# FIGURE 5 -------------------------------------------------------------------

$(OUTPUT)/reach_complexes_raw.csv: $(DATA)/bioexp_reach.pkl
	python -m bioexp.figures.figure5.validate_complex \
                             $(DATA)/bioexp_reach.pkl \
                             $(OUTPUT)/reach_complexes_raw.tsv

# KORKUT STMTS TO CHECK ------------------------------------------------------
$(OUTPUT)/bioexp_data_stmts.pkl: $(DATA)/Korkut\ et\ al.\ Data\ 05122017.xlsx
	python -u -m bioexp.explanation.make_stmts_for_checking $@

# KORKUT_PYSB ----------------------------------------------------------------
# Building the Pysb Model
$(OUTPUT)/bioexp_test_before_pysb.pkl: $(OUTPUT)/bioexp_signor.pkl
	python -u -m bioexp.explanation.assemble_pysb preprocess_stmts \
        signor test_before_pysb

$(OUTPUT)/bioexp_test_pysb_model.pkl: $(OUTPUT)/bioexp_test_before_pysb.pkl
	python -u -m bioexp.explanation.assemble_pysb assemble_pysb \
        test_before_pysb test_pysb_model

# Results
$(OUTPUT)/bioexp_test_paths.pkl: \
        $(OUTPUT)/bioexp_test_pysb_model.pkl \
        $(OUTPUT)/bioexp_data_stmts.pkl
	python -u -m bioexp.explanation.check_pysb_model test_pysb_model \
        data_stmts test_paths

#reading_only: \
#    $(OUTPUT)/bioexp_reading_only_preassembled.pkl \
#	$(OUTPUT)/bioexp_reading_only_pysb_model.pkl \
#    $(OUTPUT)/bioexp_reading_only_paths.pkl

#db_only: \
#    $(OUTPUT)/bioexp_db_only_preassembled.pkl \
#    $(OUTPUT)/bioexp_db_only_pysb_model.pkl \
#    $(OUTPUT)/bioexp_db_only_paths.pkl

