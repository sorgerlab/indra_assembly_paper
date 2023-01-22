# Must match config.json
OUTPUT := output
DATA := data
NET := networks
FIG2 := bioexp/figures/figure2
FIG4 := bioexp/figures/figure4
DEPLOY := ~/Dropbox/DARPA\ projects/papers/INDRA\ paper\ 2/figures/figure_panels

all: fig2 fig4

depmap: $(OUTPUT)/bioexp_signor_indranet_explainer.pkl

sample: $(OUTPUT)/bioexp_reach_sample_uncurated.pkl

deploy:
	rsync -av $(OUTPUT)/fig*.pdf $(DEPLOY)
	rsync -av $(OUTPUT)/fig*.png $(DEPLOY)

clean:
	cd $(OUTPUT); rm -rf *

fig2: $(OUTPUT)/fig2_evidence_distribution.pdf \
      $(OUTPUT)/fig2_stmt_counts_before_pa.pdf

fig4: \
    $(OUTPUT)/fig4_reach_curve.pdf \
    $(OUTPUT)/fig4_reach_model_fits.pdf \
    $(OUTPUT)/fig4_sparser_model_fits.pdf \
    $(OUTPUT)/fig4_medscan_model_fits.pdf \
    $(DATA)/curation/multireader_curation_dataset.pkl \

belief_fitting: $(OUTPUT)/bioexp_multi_src_results.pkl

# Change this later to point to the right stmt pickles
depmap: $(OUTPUT)/bioexp_signor_indranet.pkl

# MAKEFILE GRAPH
graph: makegraph.pdf

makegraph.pdf: makegraph.dot
	dot -T pdf makegraph.dot -o makegraph.pdf

makegraph.dot: Makefile
	make -Bnd | make2graph > makegraph.dot



# STMT SAMPLE FOR CURATION --------------------------------------------------
$(OUTPUT)/bioexp_reach_sample_uncurated.pkl: $(DATA)/bioexp_asmb_preassembled.pkl
	python -m bioexp.curation.sample $< 200 1 10 $(OUTPUT) \
                             reach sparser rlimsp isi medscan trips

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

# Compiled curation dataset for training sklearn models
$(DATA)/curation/multireader_curation_dataset.pkl: $(DATA)/bioexp_asmb_preassembled.pkl
	python -m bioexp.curation.group_curations $(DATA)/curation

# DEPMAP ----------------------------------------------------------------------
$(OUTPUT)/bioexp_signor_indranet.pkl: $(OUTPUT)/bioexp_signor.pkl
	python -m bioexp.depmap.stmts_to_indranet signor


$(OUTPUT)/bioexp_explainer.pkl: \
                              $(OUTPUT)/bioexp_signor_indranet.pkl
	python -m bioexp.depmap.get_explanations signor_indranet combined_z_score
