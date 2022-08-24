# Must match config.json
OUTPUT := output
DATA := data
NET := networks
FIG1 := bioexp/figures/figure1
FIG2 := bioexp/figures/figure2
FIG4 := bioexp/figures/figure4
DEPLOY := ~/Dropbox/DARPA\ projects/papers/INDRA\ paper\ 2/figures/figure_panels

all: fig2 fig4 fig5 korkut_pysb

depmap: $(OUTPUT)/bioexp_signor_indranet_explainer.pkl

sample: $(OUTPUT)/bioexp_reach_sample_uncurated.pkl

deploy:
	rsync -av $(OUTPUT)/fig*.pdf $(DEPLOY)
	rsync -av $(OUTPUT)/fig*.png $(DEPLOY)

clean:
	cd $(OUTPUT); rm -rf *

#fig1: $(OUTPUT)/fig1_pc_egfr_mapk1_paths.txt

fig2: $(OUTPUT)/fig2_evidence_distribution.pdf \
      $(OUTPUT)/fig2_stmt_counts_before_pa.pdf

fig4: \
    $(OUTPUT)/fig4_reach_curve.pdf \
    $(OUTPUT)/fig4_belief_surface.pdf \
    $(OUTPUT)/fig4_reach_model_fits.pdf \
    $(OUTPUT)/fig4_sparser_model_fits.pdf \
    $(OUTPUT)/fig4_medscan_model_fits.pdf \
    $(OUTPUT)/curation_dataset.pkl \
    $(OUTPUT)/fig4_ipynb_overlap_upset_linear

fig5: $(OUTPUT)/reach_complexes_raw.tsv

korkut_pysb: $(OUTPUT)/bioexp_db_only_paths.pkl

belief_fitting: $(OUTPUT)/bioexp_multi_src_results.pkl

# Change this later to point to the right stmt pickles
depmap: $(OUTPUT)/bioexp_signor_indranet.pkl

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
	python -m bioexp.figures.figure1.find_paths parse_pc \
        $(DATA)/PathwayCommons9.All.hgnc.txt \
        $(OUTPUT)/prior_genes.txt \
        $@

# STMT SAMPLE FOR CURATION --------------------------------------------------
$(OUTPUT)/bioexp_reach_sample_uncurated.pkl: $(DATA)/bioexp_asmb_preassembled.pkl
	python -m bioexp.curation.sample $< 200 1 10 $(OUTPUT) \
                             reach sparser rlimsp isi medscan trips

# FIGURE 1 -------------------------------------------------------------------

# Example pathfinding output over Pathway Commons
$(OUTPUT)/fig1_pc_egfr_mapk1_paths.txt: $(OUTPUT)/pc_multidigraph.pkl \
                                       $(FIG1)/find_paths.py
	python $(FIG1)/find_paths.py find_paths 

# FIGURE 2 -------------------------------------------------------------------

$(OUTPUT)/fig2_evidence_distribution.pdf: \
        $(DATA)/bioexp_asmb_preassembled.pkl \
        $(FIG2)/preassembly_stats.py
	python -m bioexp.figures.figure2.preassembly_stats

$(OUTPUT)/fig2_stmt_counts_before_pa.pdf: \
        $(DATA)/bioexp_asmb_preassembled.pkl \
        $(OUTPUT)/fig2_stmt_counts.txt \
        $(FIG2)/plot_stmt_counts.py
	python -m bioexp.figures.figure2.plot_stmt_counts

# FIGURE 4 -------------------------------------------------------------------
# Evidence distributions for the fits
$(OUTPUT)/bioexp_%_stmt_evidence_distribution.json:
	python -u -m bioexp.curation.get_ev_distro $*

# Run model fits (REACH)
$(OUTPUT)/fig4_model_fit_results_reach.pkl: \
    $(DATA)/curation/bioexp_reach_sample_tsv.pkl \
    $(DATA)/curation/bioexp_reach_sample_uncurated_19-12-14.pkl \
    $(DATA)/curation/bioexp_reach_sample_uncurated_20-02-19.pkl \
    $(OUTPUT)/bioexp_reach_stmt_evidence_distribution.json
	python -u -m bioexp.curation.process_curations reach $(OUTPUT)

# Run model fits (other readers)
$(OUTPUT)/fig4_model_fit_results_%.pkl: \
    $(DATA)/curation/bioexp_%_sample_uncurated.pkl \
    $(OUTPUT)/bioexp_%_stmt_evidence_distribution.json
	python -u -m bioexp.curation.process_curations $* $(OUTPUT)

# Correctness curves for each reader (no model fits)
$(OUTPUT)/fig4_reach_curve.pdf: $(FIG4)/curated_correctness.py
	python -m bioexp.figures.figure4.curated_correctness $(OUTPUT)

# Model fit plots
$(OUTPUT)/fig4_%_model_fits.pdf: \
        $(OUTPUT)/fig4_model_fit_results_%.pkl
	python -m bioexp.figures.figure4.model_fit_plots $< $* $(OUTPUT)

# Belief surface
$(OUTPUT)/fig4_belief_surface.pdf: $(FIG4)/belief_surface.py
	python -m bioexp.figures.figure4.belief_surface

# Compiled curation dataset for training sklearn models
$(OUTPUT)/curation_dataset.pkl: $(DATA)/bioexp_asmb_preassembled.pkl
	python -m bioexp.curation.group_curations $(OUTPUT)

# Various other downstreams in the notebook
# $(OUTPUT)/curation_dataset.pkl \
# $(DATA)/bioexp_asmb_preassembled.pkl
# (Check other dependencies in the notebook)

$(OUTPUT)/fig4_ipynb_overlap_upset_linear: \
        $(OUTPUT)/curation_dataset.pkl $(DATA)/bioexp_asmb_preassembled.pkl
	jupyter nbconvert --to notebook --ExecutePreprocessor.timeout=180 \
        --execute notebooks/Reader_overlap_and_error_analysis.ipynb


# FIGURE 5 -------------------------------------------------------------------

$(OUTPUT)/reach_complexes_raw.tsv: $(DATA)/bioexp_reach.pkl
	python -m bioexp.figures.figure5.validate_complex \
                             $(DATA)/bioexp_reach.pkl \
                             $(OUTPUT)/reach_complexes_raw.tsv

# KORKUT STMTS TO CHECK ------------------------------------------------------
$(OUTPUT)/bioexp_data_stmts.pkl: $(DATA)/Korkut\ et\ al.\ Data\ 05122017.xlsx
	python -u -m bioexp.explanation.make_stmts_for_checking $@

# KORKUT_PYSB ----------------------------------------------------------------
# Building the Pysb Model
$(OUTPUT)/bioexp_db_only_before_pysb.pkl: \
        $(OUTPUT)/bioexp_db_only_reduce_mods.pkl \
        $(OUTPUT)/prior_genes.txt \
        $(DATA)/antibody_site_map.csv
	python -u -m bioexp.explanation.assemble_pysb preprocess_stmts \
        db_only_reduce_mods \
        db_only_before_pysb \
        $(OUTPUT)/prior_genes.txt \
        $(DATA)/antibody_site_map.csv

$(OUTPUT)/bioexp_db_only_pysb_model.pkl: \
        $(OUTPUT)/bioexp_db_only_before_pysb.pkl \
        $(DATA)/Korkut\ et\ al.\ Data\ 05122017.xlsx \
        $(DATA)/antibody_site_map.csv
	python -u -m bioexp.explanation.assemble_pysb assemble_pysb \
        db_only_before_pysb \
        db_only_pysb_model \
        $(DATA)/Korkut\ et\ al.\ Data\ 05122017.xlsx \
        $(DATA)/antibody_site_map.csv

# Results
$(OUTPUT)/bioexp_db_only_paths.pkl: \
        $(OUTPUT)/bioexp_db_only_pysb_model.pkl \
        $(OUTPUT)/bioexp_data_stmts.pkl \
        $(DATA)/Korkut\ et\ al.\ Data\ 05122017.xlsx \
        $(DATA)/antibody_site_map.csv
	python -u -m bioexp.explanation.check_pysb_model \
        db_only_pysb_model \
        data_stmts \
        db_only_paths \
        $(DATA)/Korkut\ et\ al.\ Data\ 05122017.xlsx \
        $(DATA)/antibody_site_map.csv

#reading_only: \
#    $(OUTPUT)/bioexp_reading_only_preassembled.pkl \
#	$(OUTPUT)/bioexp_reading_only_pysb_model.pkl \
#    $(OUTPUT)/bioexp_reading_only_paths.pkl

#db_only: \
#    $(OUTPUT)/bioexp_db_only_preassembled.pkl \
#    $(OUTPUT)/bioexp_db_only_pysb_model.pkl \
#    $(OUTPUT)/bioexp_db_only_paths.pkl

# DEPMAP ----------------------------------------------------------------------
$(OUTPUT)/bioexp_signor_indranet.pkl: $(OUTPUT)/bioexp_signor.pkl
	python -m bioexp.depmap.stmts_to_indranet signor


$(OUTPUT)/bioexp_explainer.pkl: \
                              $(OUTPUT)/bioexp_signor_indranet.pkl
	python -m bioexp.depmap.get_explanations signor_indranet combined_z_score
