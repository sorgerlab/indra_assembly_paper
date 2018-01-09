BUILD := build
DATA := data

clean:
	cd $(BUILD); rm -rf *

$(BUILD)/prior_genes.txt: process_data.py \
                          $(DATA)/Korkut\ et\ al.\ Data\ 05122017.xlsx \
                          $(DATA)/ras_pathway_proteins.csv \
                          $(DATA)/drug_grounding.csv
	python process_data.py

