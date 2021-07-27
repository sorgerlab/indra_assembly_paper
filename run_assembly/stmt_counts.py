from util import pklload
from indra.tools import assemble_corpus as ac
from collections import Counter

# Reader sources
readers = ['reach', 'sparser', 'medscan', 'rlimsp',
           'isi', 'trips']
# DB sources
dbs = ['pathway_commons', 'hprd', 'bel', 'signor',
       'trrust', 'cbn', 'phosphosite']
# All sources
sources = readers + dbs

# Assembly intermediates
asmb = ['all_raw', 'filter_no_hypothesis', 'map_grounding',
        'filter_grounded_only', 'filter_genes_only', 'filter_human_only',
        'asmb_map_sequence', 'asmb_preassembled']
# All together now
filenames = sources + asmb

if __name__ == '__main__':
    stmt_counts = []
    for filename in filenames:
        # Load the pickle file
        print("Loading %s" % filename)
        if filename in ('reach', 'sparser', 'pathway_commons'):
            stmts = ac.load_statements('../data/bioexp_%s.pkl' % filename)
        else:
            stmts = pklload(filename)
        # We get the CBN UUIDs here so we can figure out later which
        # evidences are from CBN later
        if filename == 'cbn':
            cbn_uuids = {stmt.uuid for stmt in stmts}
        # Tabulate evidence source combinations
        ev_types = []
        for s in stmts:
            source_apis = []
            for e in s.evidence:
                # We use the evidence source_api by default
                source_api = e.source_api
                # We have to do some special handling for CBN here, if the
                # statements are already assembled, we have to look at
                # evidence prior UUID to check if they are originally from CBN
                if filename == 'asmb_preassembled':
                    if set(e.annotations.get('prior_uuids', set())) & cbn_uuids:
                        source_api = 'cbn'
                # Otherwise we look at the statement's UUID (only available
                # here before assembly) to see if it's originally from CBN
                else:
                    if s.uuid in cbn_uuids:
                        source_api = 'cbn'
                # We can now check if the source sub ID is phosphositeplus
                # which we need to separate from biopax
                if e.annotations.get('source_sub_id') == 'phosphositeplus':
                    source_api = 'phosphosite'
                source_apis.append(source_api)
            source_apis = frozenset(source_apis)
            ev_types.append(source_apis)
        ev_ctr = Counter(ev_types)
        print("%s stmts in %s" % (len(stmts), filename))
        stmt_counts.append((filename, len(stmts), ev_ctr))
        print(stmt_counts)
        del stmts

    # Write to TSV file with comma-separated sources in each row
    with open('../output/fig2_stmt_counts.txt', 'wt') as f:
        for filename, count, ev_ctr in stmt_counts:
            f.write('%s\ttotal\t%s\n' % (filename, count))
            for ev_set, ev_count in ev_ctr.items():
                ev_key = ','.join(ev_set)
                f.write('%s\t%s\t%s\n' % (filename, ev_key, ev_count))

