"""The purpose of this script is to identify papers that are both relevant
to the explaantion task, and are enriched in mechanisms, as highlighted by
reading with multiple reading systems. This information can be used to
prioritize what content to read with TRIPS."""

import os
import sys
import glob
import boto3
sys.path.append('..')
sys.path.append('../run_assembly')
import process_data
import util
import tqdm
import pickle
from indra.tools import assemble_corpus as ac
from collections import Counter
from functools import lru_cache

s3 = boto3.client('s3')


data = process_data.read_data()
gene_names = process_data.get_all_gene_names(data)
stmts = util.pklload('asmb_preassembled')
stmts = ac.filter_gene_list(stmts, gene_names, policy='one')

@lru_cache(maxsize=1000000)
def check_abstract_only(pmid):
    res = s3.list_objects_v2(Bucket='bigmech', Prefix='papers/PMID%s/fulltext' % pmid)
    if res['KeyCount'] == 0:
        return True
    else:
        return False


# Get all PMIDs which have reader evidence
readers = {'reach', 'sparser', 'medscan'}
pmids = []
for stmt in tqdm.tqdm(stmts):
    # We take only papers from which there was a reader extraction
    for ev in stmt.evidence:
        if ev.source_api in readers:
            if ev.pmid and check_abstract_only(ev.pmid):
                pmids.append(ev.pmid)

pmid_counts = Counter(pmids)
with open('pmid_counts.pkl', 'wb') as fh:
    pickle.dump(pmid_counts, fh)
