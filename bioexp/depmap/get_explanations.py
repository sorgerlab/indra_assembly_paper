import pandas as pd
from depmap_analysis.scripts.depmap_script2 import match_correlations
from bioexp.util import pklload


corr_file = "/Users/johnbachman/Dropbox/1johndata/knowledge file/Biology/" \
            "Research/Big Mechanism/depmap_analysis/notebooks/data/19q4/" \
            "dep_z.h5"

corr = pd.read_hdf(corr_file)
inet = pklload('signor_indranet')

