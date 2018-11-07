"""The purpose of this script is to identify papers that are both relevant
to the explaantion task, and are enriched in mechanisms, as highlighted by
reading with multiple reading systems. This information can be used to
prioritize what content to read with TRIPS."""

import sys
sys.path.append('..')
sys.path.append('../run_assembly')
import process_data
import util
from indra.tools import assemble_corpus as ac
from collections import Counter


data = process_data.read_data()
gene_names = process_data.get_all_gene_names(data)
stmts = util.pklload('filter_human_only')
stmts = ac.filter_gene_list(gene_names, policy='one')

readers = ['reach', 'sparser', 'medscan']

pmids = []
for stmt in stmts:
    for ev in stmt.evidence:
        if ev.source in readers:
            pmids.append(ev.pmid)

pmid_counts = Counter(pmids)
