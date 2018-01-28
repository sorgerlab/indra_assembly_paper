OUTPUT := output
DATA := data

clean:
	cd $(BUILD); rm -rf *

all: $(OUTPUT)/bioexp_all_raw.pkl

# DATA -----------------------------------------------------------------------

$(DATA)/%:
	python -m bioexp.transfer_s3 get $@ $(DATA)

# OUTPUT --------------------------------------------------------------

# The list of prior genes from the data and related sources

$(OUTPUT)/bioexp_bel.pkl: \
        $(DATA)/large_corpus-20170611.bel \
        process_bel.py
	python process_bel.py

$(OUTPUT)/bioexp_all_raw.pkl: \
        $(OUTPUT)/bioexp_bel.pkl \
        $(DATA)/bioexp_sparser.pkl \
        $(DATA)/bioexp_biopax.pkl \
        assemble_sources.py
	python assemble_sources.py
