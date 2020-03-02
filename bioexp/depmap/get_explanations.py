import pandas as pd
from depmap_analysis.scripts.depmap_script2 import match_correlations
from bioexp.util import pklload


corr_file = "/Users/johnbachman/Dropbox/1johndata/knowledge file/Biology/" \
            "Research/Big Mechanism/depmap_analysis/notebooks/data/19q4/" \
            "dep_z.h5"

print(f"Loading {corr_file}")
corr = pd.read_hdf(corr_file)
inet = pklload('signor_indranet')
inet_dg = inet.to_digraph()
print("Sampling")
samp = corr.sample(1000, axis=0).sample(1000, axis=1)
print(f"match_correlations")
res = match_correlations(samp, inet_dg)

