OUTPUT := output
DATA := ../data
NET := ../networks
FIG1 := ../bioexp/figures/figure1
FIG2 := ../bioexp/figures/figure2

clean:
	cd $(BUILD); rm -rf *

all: $(OUTPUT)/bioexp_bel.pkl

# DATA -----------------------------------------------------------------------

$(DATA)/%:
	python -m bioexp.transfer_s3 get $@

# OUTPUT --------------------------------------------------------------

# The list of prior genes from the data and related sources

$(OUTPUT)/bioexp_bel.pkl: \
        $(DATA)/large_corpus-20170611.bel \
        process_bel.py
	python process_bel.py


