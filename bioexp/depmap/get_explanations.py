import sys
from pathlib import Path

import pandas as pd

from bioexp.util import pklload, based, basen, prefixed_file
from depmap_analysis.scripts.depmap_script2 import main

inet_file = sys.argv[1]
corr_file = prefixed_file(sys.argv[2], 'h5')

print(f"Loading {corr_file}")
corr = pd.read_hdf(corr_file)
inet = pklload(inet_file)
inet_dg = inet.to_digraph()
print("Sampling")
samp = corr.sample(3000, axis=0)
samp = samp.filter(list(samp.index), axis=1)
print(f"match_correlations")
res = main(indra_net=inet_dg,
           sd_range=(2, 3),
           outname=str(Path(based).joinpath(basen + '_explainer').absolute()),
           graph_type='unsigned',
           z_score=corr,
           pb_model=None,
           pb_node_mapping=None,
           ignore_list=None,
           info={'test': 'test get_explanations'},
           indra_date='20200400',
           depmap_date='19Q4')
