OUTPUT := output
DATA := data

all: $(OUTPUT)/bioexp_preassembled.pkl $(OUTPUT)/fig2_stmt_counts.txt

clean:
	cd $(OUTPUT); rm -rf *

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
$(OUTPUT)/bioexp_biopax_fixed.pkl: $(DATA)/bioexp_biopax.pkl
	python -u fix_biopax.py

# ASSEMBLY ------------------------------------------------------------------

# Load the statements
$(OUTPUT)/bioexp_all_raw.pkl: \
        $(OUTPUT)/bioexp_bel.pkl \
        $(OUTPUT)/bioexp_biopax_fixed.pkl \
        $(DATA)/bioexp_reach.pkl \
        $(DATA)/bioexp_sparser.pkl
	python -u assemble_sources.py load_stmts

# Sample subset (for testing--FIXME)
#$(OUTPUT)/bioexp_all_raw_sample.pkl: $(OUTPUT)/bioexp_all_raw.pkl
#	python -u assemble_sources.py sample

# Filter no hypotheses
$(OUTPUT)/bioexp_filter_no_hypothesis.pkl: $(OUTPUT)/bioexp_all_raw.pkl
	python -u assemble_sources.py filter_no_hypothesis

# Map grounding
$(OUTPUT)/bioexp_map_grounding.pkl: \
        $(OUTPUT)/bioexp_filter_no_hypothesis.pkl
	python -u assemble_sources.py map_grounding

# Filter grounded only
$(OUTPUT)/bioexp_filter_grounded_only.pkl: \
        $(OUTPUT)/bioexp_map_grounding.pkl
	python -u assemble_sources.py filter_grounded_only

# Filter genes only
$(OUTPUT)/bioexp_filter_genes_only.pkl: \
        $(OUTPUT)/bioexp_filter_grounded_only.pkl
	python -u assemble_sources.py filter_genes_only

# Filter human only
$(OUTPUT)/bioexp_filter_human_only.pkl: \
        $(OUTPUT)/bioexp_filter_genes_only.pkl
	python -u assemble_sources.py filter_human_only

# Expand families
$(OUTPUT)/bioexp_expand_families.pkl: \
        $(OUTPUT)/bioexp_filter_human_only.pkl
	python -u assemble_sources.py expand_families

# Filter to gene list
$(OUTPUT)/bioexp_filter_gene_list.pkl: \
        $(OUTPUT)/bioexp_expand_families.pkl \
        ../build/prior_genes.txt
	python -u assemble_sources.py filter_gene_list

# Map sequences
$(OUTPUT)/bioexp_map_sequence.pkl: \
        $(OUTPUT)/bioexp_filter_gene_list.pkl
	python -u assemble_sources.py map_sequence

# Run preassembly
$(OUTPUT)/bioexp_preassembled.pkl: \
        $(OUTPUT)/bioexp_map_sequence.pkl
	python -u assemble_sources.py run_preassembly

# -- ANALYSIS/STATISTICS -----------------------------------------------------

$(OUTPUT)/fig2_stmt_counts.txt: \
        $(OUTPUT)/bioexp_reach.pkl \
        $(OUTPUT)/bioexp_sparser.pkl \
        $(OUTPUT)/bioexp_bel.pkl \
        $(OUTPUT)/bioexp_fixed.pkl \
        $(OUTPUT)/bioexp_all_raw.pkl \
        $(OUTPUT)/bioexp_filter_no_hypothesis.pkl \
        $(OUTPUT)/bioexp_map_grounding.pkl \
        $(OUTPUT)/bioexp_filter_grounded_only.pkl \
        $(OUTPUT)/bioexp_filter_genes_only.pkl \
        $(OUTPUT)/bioexp_filter_human_only.pkl \
        $(OUTPUT)/bioexp_expand_families.pkl \
        $(OUTPUT)/bioexp_filter_gene_list.pkl \
        $(OUTPUT)/bioexp_map_sequence.pkl \
        $(OUTPUT)/bioexp_preassembled.pkl
	python -u stmt_counts.py



