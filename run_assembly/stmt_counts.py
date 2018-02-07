from util import pklload
from indra.tools import assemble_corpus as ac
from collections import Counter

filenames = ['reach', 'sparser', 'bel', 'signor', 'biopax_fixed', 'all_raw',
             'filter_no_hypothesis', 'map_grounding', 'filter_grounded_only',
             'filter_genes_only', 'filter_human_only', 'expand_families',
             'filter_gene_list', 'map_sequence', 'preassembled',
             'filter_belief', 'filter_top_level', 'filter_enzyme_kinase',
             'filter_mod_nokinase', 'reduce_activities', 'reduce_mods']

stmt_counts = []
for filename in filenames:
    # Load the pickle file
    print("Loading %s" % filename)
    if filename in ('reach', 'sparser', 'biopax_fixed'):
        stmts = ac.load_statements('data/bioexp_%s.pkl' % filename)
    else:
        stmts = pklload(filename)
    # Tabulate evidence source combinations
    ev_types = []
    for s in stmts:
        source_apis = []
        for e in s.evidence:
            source_apis.append(e.source_api)
        source_apis = frozenset(source_apis)
        ev_types.append(source_apis)
    ev_ctr = Counter(ev_types)
    print("%s stmts in %s" % (len(stmts), filename))
    stmt_counts.append((filename, len(stmts), ev_ctr))

# Write to TSV file with comma-separated sources in each row
with open('output/fig2_stmt_counts.txt', 'wt') as f:
    for filename, count, ev_ctr in stmt_counts:
        f.write('%s\ttotal\t%s\n' % (filename, count))
        for ev_set, ev_count in ev_ctr.items():
            ev_key = ','.join(ev_set)
            f.write('%s\t%s\t%s\n' % (filename, ev_key, ev_count))

