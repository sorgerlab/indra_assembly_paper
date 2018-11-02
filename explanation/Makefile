OUTPUT := ~/data/bioexp
DATA := ../run_assembly/data

all: assembly explanation stats

reading_only: \
    $(OUTPUT)/bioexp_reading_only_asmb_preassembled.pkl \
    $(OUTPUT)/bioexp_reach_only_preassembled.pkl \
    $(OUTPUT)/bioexp_medscan_only_preassembled.pkl \
    $(OUTPUT)/bioexp_sparser_only_preassembled.pkl

explanation: \
    $(OUTPUT)/bioexp_reduce_mods.pkl \
    $(OUTPUT)/bioexp_db_only_reduce_mods.pkl \
    $(OUTPUT)/bioexp_reading_only_reduce_mods.pkl

db_only: \
	$(OUTPUT)/bioexp_db_only_before_pysb.pkl \
	$(OUTPUT)/bioexp_db_only_pysb_model.pkl

stats: $(OUTPUT)/fig2_stmt_counts.txt

clean:
	cd $(OUTPUT); rm -rf *

s3:
	aws s3 cp --recursive output s3://bioexp-paper


# -- ANALYSIS/STATISTICS -----------------------------------------------------

$(OUTPUT)/bioexp_db_only_before_pysb.pkl: \
        $(OUTPUT)/bioexp_db_only_reduce_mods.pkl
	python -u assemble_pysb.py preprocess_statements \
        db_only_reduce_mods db_only_before_pysb

$(OUTPUT)/bioexp_db_only_pysb_model.pkl: \
        $(OUTPUT)/bioexp_db_only_before_pysb.pkl
	python -u assemble_pysb.py assemble_pysb \
        db_only_before_pysb db_only_pysb_model
