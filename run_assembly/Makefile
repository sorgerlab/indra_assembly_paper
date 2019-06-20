OUTPUT := ~/data/bioexp
DATA := data

all: assembly explanation stats

assembly: \
    $(OUTPUT)/bioexp_asmb_preassembled.pkl \
    $(OUTPUT)/bioexp_db_only_asmb_preassembled.pkl \
    reading_only_assembly

reading_only_assembly: \
    $(OUTPUT)/bioexp_reading_only_asmb_preassembled.pkl \
    $(OUTPUT)/bioexp_reach_only_preassembled.pkl \
    $(OUTPUT)/bioexp_medscan_only_preassembled.pkl \
    $(OUTPUT)/bioexp_sparser_only_preassembled.pkl \
    $(OUTPUT)/bioexp_trips_only_preassembled.pkl

explanation: \
    $(OUTPUT)/bioexp_reduce_mods.pkl \
    $(OUTPUT)/bioexp_db_only_reduce_mods.pkl \
    $(OUTPUT)/bioexp_reading_only_reduce_mods.pkl

db_only: \
    $(OUTPUT)/bioexp_db_only_asmb_preassembled.pkl \
    $(OUTPUT)/bioexp_db_only_reduce_mods.pkl

stats: $(OUTPUT)/fig2_stmt_counts.txt

clean:
	cd $(OUTPUT); rm -rf *

s3:
	aws s3 cp --recursive output s3://bioexp-paper

# MAKEFILE GRAPH
graph: makegraph.pdf

makegraph.pdf: makegraph.dot
	dot -T pdf makegraph.dot -o makegraph.pdf

makegraph.dot: Makefile
	make -Bnd | make2graph > makegraph.dot


# DATA -----------------------------------------------------------------------

$(DATA)/%:
	python -u -m bioexp.transfer_s3 get $@ $(DATA)

# PREPROCESSING --------------------------------------------------------------

# Process the BEL large corpus
$(OUTPUT)/bioexp_bel.pkl: \
        $(DATA)/large_corpus-20170611.bel \
        process_bel.py
	python -u process_bel.py

# Fix the biopax pickle
#$(OUTPUT)/bioexp_biopax_fixed.pkl: $(DATA)/bioexp_biopax.pkl
#	python -u fix_biopax.py

# Download Signor statements
$(OUTPUT)/bioexp_signor.pkl:
	python -u assemble_sources.py signor

# Download Biogrid statements
$(OUTPUT)/bioexp_biogrid.pkl:
	python -u get_biogrid_data.py

# FULL ASSEMBLY PIPELINE (DBs + Reading) -------------------------------------

# Load the statements
$(OUTPUT)/bioexp_all_raw.pkl: \
        $(OUTPUT)/bioexp_signor.pkl \
        $(OUTPUT)/bioexp_bel.pkl \
        $(OUTPUT)/bioexp_biopax.pkl \
        $(OUTPUT)/bioexp_hprd.pkl \
        $(OUTPUT)/bioexp_trrust.pkl \
        $(OUTPUT)/bioexp_cbn.pkl \
        $(OUTPUT)/bioexp_phosphosite.pkl \
        $(OUTPUT)/bioexp_reach.pkl \
        $(OUTPUT)/bioexp_sparser.pkl \
        $(OUTPUT)/bioexp_trips.pkl \
        $(OUTPUT)/bioexp_rlimsp.pkl \
        $(OUTPUT)/bioexp_medscan.pkl \
        $(OUTPUT)/bioexp_isi.pkl
	python -u assemble_sources.py load_stmts all_raw \
        $(OUTPUT)/bioexp_signor.pkl \
        $(OUTPUT)/bioexp_bel.pkl \
        $(OUTPUT)/bioexp_biopax.pkl \
        $(OUTPUT)/bioexp_hprd.pkl \
        $(OUTPUT)/bioexp_trrust.pkl \
        $(OUTPUT)/bioexp_cbn.pkl \
        $(OUTPUT)/bioexp_phosphosite.pkl \
        $(OUTPUT)/bioexp_reach.pkl \
        $(OUTPUT)/bioexp_sparser.pkl \
        $(OUTPUT)/bioexp_trips.pkl \
        $(OUTPUT)/bioexp_rlimsp.pkl \
        $(OUTPUT)/bioexp_medscan.pkl \
        $(OUTPUT)/bioexp_isi.pkl

# Filter no hypotheses
$(OUTPUT)/bioexp_filter_no_hypothesis.pkl: $(OUTPUT)/bioexp_all_raw.pkl
	python -u assemble_sources.py filter_no_hypothesis \
        all_raw filter_no_hypothesis

# Map grounding
$(OUTPUT)/bioexp_map_grounding.pkl: $(OUTPUT)/bioexp_filter_no_hypothesis.pkl
	python -u assemble_sources.py map_grounding \
        filter_no_hypothesis map_grounding

# Filter grounded only
$(OUTPUT)/bioexp_filter_grounded_only.pkl: $(OUTPUT)/bioexp_map_grounding.pkl
	python -u assemble_sources.py filter_grounded_only \
        map_grounding filter_grounded_only

# Filter genes only
$(OUTPUT)/bioexp_filter_genes_only.pkl: \
        $(OUTPUT)/bioexp_filter_grounded_only.pkl
	python -u assemble_sources.py filter_genes_only \
        filter_grounded_only filter_genes_only

# Filter human only
$(OUTPUT)/bioexp_filter_human_only.pkl: $(OUTPUT)/bioexp_filter_genes_only.pkl
	python -u assemble_sources.py filter_human_only \
        filter_genes_only filter_human_only

# ---------------------------
# Branch -> Assembly only
# Map sequences
$(OUTPUT)/bioexp_asmb_map_sequence.pkl: \
        $(OUTPUT)/bioexp_filter_human_only.pkl
	python -u assemble_sources.py map_sequence \
        filter_human_only asmb_map_sequence

# Run preassembly
$(OUTPUT)/bioexp_asmb_preassembled.pkl: \
        $(OUTPUT)/bioexp_asmb_map_sequence.pkl
	python -u assemble_sources.py preassembled \
        asmb_map_sequence asmb_preassembled

# ---------------------------
# Branch -> Explanation/modeling (expand families)
# Expand families
$(OUTPUT)/bioexp_expand_families.pkl: $(OUTPUT)/bioexp_filter_human_only.pkl
	python -u assemble_sources.py expand_families \
        filter_human_only expand_families

# Filter to gene list
$(OUTPUT)/bioexp_filter_gene_list.pkl: \
        $(OUTPUT)/bioexp_expand_families.pkl \
        ../build/prior_genes.txt
	python -u assemble_sources.py filter_gene_list \
        expand_families filter_gene_list

# Map sequences
$(OUTPUT)/bioexp_map_sequence.pkl: $(OUTPUT)/bioexp_filter_gene_list.pkl
	python -u assemble_sources.py map_sequence filter_gene_list map_sequence

# Run preassembly
$(OUTPUT)/bioexp_preassembled.pkl: $(OUTPUT)/bioexp_map_sequence.pkl
	python -u assemble_sources.py preassembled map_sequence preassembled

# Filter out high belief statements
$(OUTPUT)/bioexp_filter_belief.pkl: $(OUTPUT)/bioexp_preassembled.pkl
	python -u assemble_sources.py filter_belief preassembled filter_belief

# Filter out all stmts except top-level
$(OUTPUT)/bioexp_filter_top_level.pkl: $(OUTPUT)/bioexp_filter_belief.pkl
	python -u assemble_sources.py filter_top_level \
        filter_belief filter_top_level

# Filter phosphorylations to stmts with known kinases as enzymes
$(OUTPUT)/bioexp_filter_enzyme_kinase.pkl: $(OUTPUT)/bioexp_filter_top_level.pkl
	python -u assemble_sources.py filter_enzyme_kinase \
        filter_top_level filter_enzyme_kinase

# Filter out non-phos mods with kinases as enzymes
$(OUTPUT)/bioexp_filter_mod_nokinase.pkl: \
        $(OUTPUT)/bioexp_filter_enzyme_kinase.pkl
	python -u assemble_sources.py filter_mod_nokinase \
        filter_enzyme_kinase filter_mod_nokinase

# Filter out non-transcription factors transcribing proteins
$(OUTPUT)/bioexp_filter_transcription_factor.pkl: \
        $(OUTPUT)/bioexp_filter_mod_nokinase.pkl
	python -u assemble_sources.py filter_transcription_factor \
        filter_mod_nokinase filter_transcription_factor

# Infer out redundant activities
$(OUTPUT)/bioexp_reduce_activities.pkl: \
        $(OUTPUT)/bioexp_filter_transcription_factor.pkl
	python -u assemble_sources.py reduce_activities \
        filter_transcription_factor reduce_activities

# Infer out redundant mods
$(OUTPUT)/bioexp_reduce_mods.pkl: $(OUTPUT)/bioexp_reduce_activities.pkl
	python -u assemble_sources.py reduce_mods \
        reduce_activities reduce_mods


# READING-ONLY ASSEMBLY PIPELINE (for evaluation) ----------------------------

$(OUTPUT)/bioexp_reading_only_all_raw.pkl: \
        $(OUTPUT)/bioexp_reach.pkl \
        $(OUTPUT)/bioexp_sparser.pkl \
        $(OUTPUT)/bioexp_trips.pkl \
        $(OUTPUT)/bioexp_rlimsp.pkl \
        $(OUTPUT)/bioexp_medscan.pkl \
        $(OUTPUT)/bioexp_isi.pkl

	python -u assemble_sources.py load_stmts reading_only_all_raw \
        $(OUTPUT)/bioexp_reach.pkl \
        $(OUTPUT)/bioexp_sparser.pkl \
        $(OUTPUT)/bioexp_trips.pkl \
        $(OUTPUT)/bioexp_rlimsp.pkl \
        $(OUTPUT)/bioexp_medscan.pkl \
        $(OUTPUT)/bioexp_isi.pkl

# Filter no hypotheses
$(OUTPUT)/bioexp_reading_only_filter_no_hypothesis.pkl: \
        $(OUTPUT)/bioexp_reading_only_all_raw.pkl
	python -u assemble_sources.py filter_no_hypothesis \
        reading_only_all_raw reading_only_filter_no_hypothesis

# Map grounding
$(OUTPUT)/bioexp_reading_only_map_grounding.pkl: \
        $(OUTPUT)/bioexp_reading_only_filter_no_hypothesis.pkl
	python -u assemble_sources.py map_grounding \
        reading_only_filter_no_hypothesis reading_only_map_grounding

# Filter grounded only
$(OUTPUT)/bioexp_reading_only_filter_grounded_only.pkl: \
        $(OUTPUT)/bioexp_reading_only_map_grounding.pkl
	python -u assemble_sources.py filter_grounded_only \
        reading_only_map_grounding reading_only_filter_grounded_only

# Filter genes only
$(OUTPUT)/bioexp_reading_only_filter_genes_only.pkl: \
        $(OUTPUT)/bioexp_reading_only_filter_grounded_only.pkl
	python -u assemble_sources.py filter_genes_only \
        reading_only_filter_grounded_only reading_only_filter_genes_only

# Filter human only
$(OUTPUT)/bioexp_reading_only_filter_human_only.pkl: \
        $(OUTPUT)/bioexp_reading_only_filter_genes_only.pkl
	python -u assemble_sources.py filter_human_only \
        reading_only_filter_genes_only reading_only_filter_human_only

# ---------------------------
# Branch -> Assembly only
# Map sequences
$(OUTPUT)/bioexp_reading_only_asmb_map_sequence.pkl: \
        $(OUTPUT)/bioexp_reading_only_filter_human_only.pkl
	python -u assemble_sources.py map_sequence \
        reading_only_filter_human_only reading_only_asmb_map_sequence

# Run preassembly
$(OUTPUT)/bioexp_reading_only_asmb_preassembled.pkl: \
        $(OUTPUT)/bioexp_reading_only_asmb_map_sequence.pkl
	python -u assemble_sources.py preassembled \
        reading_only_asmb_map_sequence reading_only_asmb_preassembled

# ---------------------------
# Branch -> Explanation/modeling (expand families)
# Expand families
$(OUTPUT)/bioexp_reading_only_expand_families.pkl: \
        $(OUTPUT)/bioexp_reading_only_filter_human_only.pkl
	python -u assemble_sources.py expand_families \
        reading_only_filter_human_only reading_only_expand_families

# Map sequences
$(OUTPUT)/bioexp_reading_only_map_sequence.pkl: \
        $(OUTPUT)/bioexp_reading_only_expand_families.pkl
	python -u assemble_sources.py map_sequence \
        reading_only_expand_families reading_only_map_sequence

# Run preassembly
$(OUTPUT)/bioexp_reading_only_preassembled.pkl: \
        $(OUTPUT)/bioexp_reading_only_map_sequence.pkl
	python -u assemble_sources.py preassembled \
        reading_only_map_sequence reading_only_preassembled

# Filter out high belief statements
$(OUTPUT)/bioexp_reading_only_filter_belief.pkl: \
        $(OUTPUT)/bioexp_reading_only_preassembled.pkl
	python -u assemble_sources.py filter_belief \
        reading_only_preassembled reading_only_filter_belief

# Filter out all stmts except top-level
$(OUTPUT)/bioexp_reading_only_filter_top_level.pkl: \
        $(OUTPUT)/bioexp_reading_only_filter_belief.pkl
	python -u assemble_sources.py filter_top_level \
        reading_only_filter_belief reading_only_filter_top_level

# Filter phosphorylations to stmts with known kinases as enzymes
$(OUTPUT)/bioexp_reading_only_filter_enzyme_kinase.pkl: \
        $(OUTPUT)/bioexp_reading_only_filter_top_level.pkl
	python -u assemble_sources.py filter_enzyme_kinase \
        reading_only_filter_top_level reading_only_filter_enzyme_kinase

# Filter out non-phos mods with kinases as enzymes
$(OUTPUT)/bioexp_reading_only_filter_mod_nokinase.pkl: \
        $(OUTPUT)/bioexp_reading_only_filter_enzyme_kinase.pkl
	python -u assemble_sources.py filter_mod_nokinase \
        reading_only_filter_enzyme_kinase reading_only_filter_mod_nokinase

# Filter out non-transcription factors transcribing proteins
$(OUTPUT)/bioexp_reading_only_filter_transcription_factor.pkl: \
        $(OUTPUT)/bioexp_reading_only_filter_mod_nokinase.pkl
	python -u assemble_sources.py filter_transcription_factor \
       reading_only_filter_mod_nokinase reading_only_filter_transcription_factor

# Infer out redundant activities
$(OUTPUT)/bioexp_reading_only_reduce_activities.pkl: \
        $(OUTPUT)/bioexp_reading_only_filter_transcription_factor.pkl
	python -u assemble_sources.py reduce_activities \
        reading_only_filter_transcription_factor reading_only_reduce_activities

# Infer out redundant mods
$(OUTPUT)/bioexp_reading_only_reduce_mods.pkl: \
        $(OUTPUT)/bioexp_reading_only_reduce_activities.pkl
	python -u assemble_sources.py reduce_mods \
        reading_only_reduce_activities reading_only_reduce_mods

# SOURCE SPECIFIC ASSEMBLY PIPELINES (for comparison to all reading) ----------

# REACH only
$(OUTPUT)/bioexp_reach_only.pkl: \
        $(OUTPUT)/bioexp_reading_only_asmb_map_sequence.pkl
	python -u assemble_sources.py filter_source reading_only_asmb_map_sequence \
        reach reach_only

$(OUTPUT)/bioexp_reach_only_preassembled.pkl: \
        $(OUTPUT)/bioexp_reach_only.pkl
	python -u assemble_sources.py preassembled \
        reach_only reach_only_preassembled

# MEDSCAN only
$(OUTPUT)/bioexp_medscan_only.pkl: \
        $(OUTPUT)/bioexp_reading_only_asmb_map_sequence.pkl
	python -u assemble_sources.py filter_source reading_only_asmb_map_sequence \
        medscan medscan_only

$(OUTPUT)/bioexp_medscan_only_preassembled.pkl: \
        $(OUTPUT)/bioexp_medscan_only.pkl
	python -u assemble_sources.py preassembled \
        medscan_only medscan_only_preassembled

# SPARSER only
$(OUTPUT)/bioexp_sparser_only.pkl: \
        $(OUTPUT)/bioexp_reading_only_asmb_map_sequence.pkl
	python -u assemble_sources.py filter_source reading_only_asmb_map_sequence \
        sparser sparser_only

$(OUTPUT)/bioexp_sparser_only_preassembled.pkl: \
        $(OUTPUT)/bioexp_sparser_only.pkl
	python -u assemble_sources.py preassembled \
        sparser_only sparser_only_preassembled


# TRIPS only
$(OUTPUT)/bioexp_trips_only.pkl: \
        $(OUTPUT)/bioexp_reading_only_asmb_map_sequence.pkl
	python -u assemble_sources.py filter_source reading_only_asmb_map_sequence \
        trips trips_only

$(OUTPUT)/bioexp_trips_only_preassembled.pkl: \
        $(OUTPUT)/bioexp_trips_only.pkl
	python -u assemble_sources.py preassembled \
        trips_only trips_only_preassembled

# DATABASE-ONLY ASSEMBLY PIPELINE (for evaluation) ----------------------------

$(OUTPUT)/bioexp_db_only_all_raw.pkl: \
        $(OUTPUT)/bioexp_signor.pkl \
        $(OUTPUT)/bioexp_bel.pkl \
        $(OUTPUT)/bioexp_biopax.pkl \
        $(OUTPUT)/bioexp_hprd.pkl \
        $(OUTPUT)/bioexp_trrust.pkl \
        $(OUTPUT)/bioexp_cbn.pkl \
        $(OUTPUT)/bioexp_phosphosite.pkl
	python -u assemble_sources.py load_stmts db_only_all_raw \
        $(OUTPUT)/bioexp_signor.pkl \
        $(OUTPUT)/bioexp_bel.pkl \
        $(OUTPUT)/bioexp_biopax.pkl \
        $(OUTPUT)/bioexp_hprd.pkl \
        $(OUTPUT)/bioexp_trrust.pkl \
        $(OUTPUT)/bioexp_cbn.pkl \
        $(OUTPUT)/bioexp_phosphosite.pkl

# Filter no hypotheses
$(OUTPUT)/bioexp_db_only_filter_no_hypothesis.pkl: \
        $(OUTPUT)/bioexp_db_only_all_raw.pkl
	python -u assemble_sources.py filter_no_hypothesis \
        db_only_all_raw db_only_filter_no_hypothesis

# Map grounding
$(OUTPUT)/bioexp_db_only_map_grounding.pkl: \
        $(OUTPUT)/bioexp_db_only_filter_no_hypothesis.pkl
	python -u assemble_sources.py map_grounding \
        db_only_filter_no_hypothesis db_only_map_grounding

# Filter grounded only
$(OUTPUT)/bioexp_db_only_filter_grounded_only.pkl: \
        $(OUTPUT)/bioexp_db_only_map_grounding.pkl
	python -u assemble_sources.py filter_grounded_only \
        db_only_map_grounding db_only_filter_grounded_only

# Filter genes only
$(OUTPUT)/bioexp_db_only_filter_genes_only.pkl: \
        $(OUTPUT)/bioexp_db_only_filter_grounded_only.pkl
	python -u assemble_sources.py filter_genes_only \
        db_only_filter_grounded_only db_only_filter_genes_only

# Filter human only
$(OUTPUT)/bioexp_db_only_filter_human_only.pkl: \
        $(OUTPUT)/bioexp_db_only_filter_genes_only.pkl
	python -u assemble_sources.py filter_human_only \
        db_only_filter_genes_only db_only_filter_human_only

# ---------------------------
# Branch -> Assembly only
# Map sequences
$(OUTPUT)/bioexp_db_only_asmb_map_sequence.pkl: \
        $(OUTPUT)/bioexp_db_only_filter_human_only.pkl
	python -u assemble_sources.py map_sequence \
        db_only_filter_human_only db_only_asmb_map_sequence

# Run preassembly
$(OUTPUT)/bioexp_db_only_asmb_preassembled.pkl: \
        $(OUTPUT)/bioexp_db_only_asmb_map_sequence.pkl
	python -u assemble_sources.py preassembled \
        db_only_asmb_map_sequence db_only_asmb_preassembled

# ---------------------------
# Branch -> Explanation/modeling (expand families)
# Expand families
$(OUTPUT)/bioexp_db_only_expand_families.pkl: \
        $(OUTPUT)/bioexp_db_only_filter_human_only.pkl
	python -u assemble_sources.py expand_families \
        db_only_filter_human_only db_only_expand_families

# Map sequences
$(OUTPUT)/bioexp_db_only_map_sequence.pkl: \
        $(OUTPUT)/bioexp_db_only_expand_families.pkl
	python -u assemble_sources.py map_sequence \
        db_only_expand_families db_only_map_sequence

# Run preassembly
$(OUTPUT)/bioexp_db_only_preassembled.pkl: \
        $(OUTPUT)/bioexp_db_only_map_sequence.pkl
	python -u assemble_sources.py preassembled \
        db_only_map_sequence db_only_preassembled

# Filter out high belief statements
$(OUTPUT)/bioexp_db_only_filter_belief.pkl: \
        $(OUTPUT)/bioexp_db_only_preassembled.pkl
	python -u assemble_sources.py filter_belief \
        db_only_preassembled db_only_filter_belief

# Filter out all stmts except top-level
$(OUTPUT)/bioexp_db_only_filter_top_level.pkl: \
        $(OUTPUT)/bioexp_db_only_filter_belief.pkl
	python -u assemble_sources.py filter_top_level \
        db_only_filter_belief db_only_filter_top_level

# Filter phosphorylations to stmts with known kinases as enzymes
$(OUTPUT)/bioexp_db_only_filter_enzyme_kinase.pkl: \
        $(OUTPUT)/bioexp_db_only_filter_top_level.pkl
	python -u assemble_sources.py filter_enzyme_kinase \
        db_only_filter_top_level db_only_filter_enzyme_kinase

# Filter out non-phos mods with kinases as enzymes
$(OUTPUT)/bioexp_db_only_filter_mod_nokinase.pkl: \
        $(OUTPUT)/bioexp_db_only_filter_enzyme_kinase.pkl
	python -u assemble_sources.py filter_mod_nokinase \
        db_only_filter_enzyme_kinase db_only_filter_mod_nokinase

# Filter out non-transcription factors transcribing proteins
$(OUTPUT)/bioexp_db_only_filter_transcription_factor.pkl: \
        $(OUTPUT)/bioexp_db_only_filter_mod_nokinase.pkl
	python -u assemble_sources.py filter_transcription_factor \
        db_only_filter_mod_nokinase db_only_filter_transcription_factor

# Infer out redundant activities
$(OUTPUT)/bioexp_db_only_reduce_activities.pkl: \
        $(OUTPUT)/bioexp_db_only_filter_transcription_factor.pkl
	python -u assemble_sources.py reduce_activities \
        db_only_filter_transcription_factor db_only_reduce_activities

# Infer out redundant mods
$(OUTPUT)/bioexp_db_only_reduce_mods.pkl: \
        $(OUTPUT)/bioexp_db_only_reduce_activities.pkl
	python -u assemble_sources.py reduce_mods \
        db_only_reduce_activities db_only_reduce_mods


# -- ANALYSIS/STATISTICS -----------------------------------------------------

$(OUTPUT)/fig2_stmt_counts.txt: \
        $(OUTPUT)/bioexp_signor.pkl \
        $(OUTPUT)/bioexp_bel.pkl \
        $(OUTPUT)/bioexp_biopax.pkl \
        $(OUTPUT)/bioexp_hprd.pkl \
        $(OUTPUT)/bioexp_trrust.pkl \
        $(OUTPUT)/bioexp_cbn.pkl \
        $(OUTPUT)/bioexp_phosphosite.pkl \
        $(OUTPUT)/bioexp_reach.pkl \
        $(OUTPUT)/bioexp_sparser.pkl \
        $(OUTPUT)/bioexp_trips.pkl \
        $(OUTPUT)/bioexp_rlimsp.pkl \
        $(OUTPUT)/bioexp_medscan.pkl \
        $(OUTPUT)/bioexp_isi.pkl \
        $(OUTPUT)/bioexp_all_raw.pkl \
        $(OUTPUT)/bioexp_filter_no_hypothesis.pkl \
        $(OUTPUT)/bioexp_map_grounding.pkl \
        $(OUTPUT)/bioexp_filter_grounded_only.pkl \
        $(OUTPUT)/bioexp_filter_genes_only.pkl \
        $(OUTPUT)/bioexp_filter_human_only.pkl \
        $(OUTPUT)/bioexp_expand_families.pkl \
        $(OUTPUT)/bioexp_filter_gene_list.pkl \
        $(OUTPUT)/bioexp_map_sequence.pkl \
        $(OUTPUT)/bioexp_preassembled.pkl \
        $(OUTPUT)/bioexp_filter_belief.pkl \
        $(OUTPUT)/bioexp_filter_top_level.pkl \
        $(OUTPUT)/bioexp_filter_enzyme_kinase.pkl \
        $(OUTPUT)/bioexp_filter_mod_nokinase.pkl \
        $(OUTPUT)/bioexp_filter_transcription_factor.pkl \
        $(OUTPUT)/bioexp_reduce_activities.pkl \
        $(OUTPUT)/bioexp_reduce_mods.pkl
	python -u stmt_counts.py

$(OUTPUT)/stmt_overlap_reading.pdf: \
        $(OUTPUT)/bioexp_reach_only_preassembled.pkl \
        $(OUTPUT)/bioexp_medscan_only_preassembled.pkl \
        $(OUTPUT)/bioexp_sparser_only_preassembled.pkl \
        $(OUTPUT)/bioexp_trips_only_preassembled.pkl \
        $(OUTPUT)/bioexp_reading_only_asmb_preassembled.pkl
	python -u multi_reader_analysis.py


#asdfasdf: \
#        $(OUTPUT)/bioexp_biogrid.pkl

