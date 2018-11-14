OUTPUT := ~/data/bioexp
DATA := ../run_assembly/data

all: db_only reading_only

reading_only: \
    $(OUTPUT)/bioexp_reading_only_preassembled.pkl \
	$(OUTPUT)/bioexp_reading_only_pysb_model.pkl \
    $(OUTPUT)/bioexp_reading_only_paths.pkl

db_only: \
    $(OUTPUT)/bioexp_db_only_preassembled.pkl \
    $(OUTPUT)/bioexp_db_only_pysb_model.pkl \
    $(OUTPUT)/bioexp_db_only_paths.pkl

$(OUTPUT)/bioexp_db_only_before_pysb.pkl: \
        assemble_pysb.py \
        $(OUTPUT)/bioexp_db_only_reduce_mods.pkl
	python -u assemble_pysb.py preprocess_stmts \
        db_only_reduce_mods db_only_before_pysb

$(OUTPUT)/bioexp_db_only_pysb_model.pkl: \
        $(OUTPUT)/bioexp_db_only_before_pysb.pkl
	python -u assemble_pysb.py assemble_pysb \
        db_only_before_pysb db_only_pysb_model



# -- ANALYSIS/STATISTICS -----------------------------------------------------

$(OUTPUT)/bioexp_db_only_paths.pkl: \
        $(OUTPUT)/bioexp_db_only_pysb_model.pkl
	python -u check_pysb_model.py db_only_pysb_model \
        db_only_paths
