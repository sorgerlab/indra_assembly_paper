OUTPUT := output
DATA := data

all: $(OUTPUT)/bioexp_test_paths.pkl

#all: db_only reading_only

#reading_only: \
#    $(OUTPUT)/bioexp_reading_only_preassembled.pkl \
#	$(OUTPUT)/bioexp_reading_only_pysb_model.pkl \
#    $(OUTPUT)/bioexp_reading_only_paths.pkl

#db_only: \
#    $(OUTPUT)/bioexp_db_only_preassembled.pkl \
#    $(OUTPUT)/bioexp_db_only_pysb_model.pkl \
#    $(OUTPUT)/bioexp_db_only_paths.pkl

$(OUTPUT)/bioexp_test_before_pysb.pkl: $(OUTPUT)/bioexp_signor.pkl
	python -u -m bioexp.explanation.assemble_pysb preprocess_stmts \
        signor test_before_pysb

$(OUTPUT)/bioexp_test_pysb_model.pkl: $(OUTPUT)/bioexp_test_before_pysb.pkl
	python -u -m bioexp.explanation.assemble_pysb assemble_pysb \
        test_before_pysb test_pysb_model



# -- ANALYSIS/STATISTICS -----------------------------------------------------

$(OUTPUT)/bioexp_test_paths.pkl: \
        $(OUTPUT)/bioexp_test_pysb_model.pkl
	python -u -m bioexp.explanation.check_pysb_model test_pysb_model \
        test_paths
