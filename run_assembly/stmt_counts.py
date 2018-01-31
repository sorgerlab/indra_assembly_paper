from util import pklload

filenames = ['reach', 'sparser', 'bel', 'biopax_fixed', 'all_raw',
             'filter_no_hypothesis', 'map_grounding', 'filter_grounded_only',
             'filter_genes_only', 'filter_human_only', 'expand_families',
             'filter_gene_list', 'map_sequence', 'preassembled']

stmt_counts = []
for filename in filenames:
    print("Loading %s" % filename)
    if filename in ('reach', 'sparser'):
        with open('data/bioexp_%s.pkl' % filename, 'rb') as f:
            stmts = pickle.load(f)
    else:
        stmts = pklload(filename)
    print("%s stmts in %s" % (len(stmts), filename))
    stmt_counts.append((filename, len(stmts)))

with open('fig2_stmt_counts.txt', 'wt') as f:
    for filename, count in stmt_counts:
        f.write('%s\t%s\n' % (filename, count))

