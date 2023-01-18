# This script replaces pickles of statements sampled from the assembled
# statements for curation purposes with lists of hashes of those statements.
# This allows for avoiding the use of non-public and large pickle files
# like bioexp_rlimsp_sample_uncurated.pkl and using instead just the
# hashes of the statements plus filtering of the single preassembled
# pickle.

import copy
import json
from os.path import join, dirname, abspath
from indra.tools import assemble_corpus as ac
from bioexp.curation.process_curations import reader_input, \
    load_curated_pkl_files, curation_data

print('Loading preassembled statements')
asmb_pkl = join(dirname(abspath(__file__)), '..', '..', 'data',
                'bioexp_asmb_preassembled.pkl')
all_stmts = ac.load_statements(asmb_pkl)
all_stmts_by_hash = {stmt.get_hash(): stmt for stmt in all_stmts}

for reader, data in reader_input.items():
    print('Loading %s' % reader)
    for pkl_fname in data['pkl_list']:
        print('Loading %s' % pkl_fname)
        stmts = load_curated_pkl_files([pkl_fname])
        stmt_hashes = [stmt.get_hash() for stmt in stmts]
        fname = join(curation_data, pkl_fname.replace('.pkl', '_hashes.json'))
        # Now write the hashes into a JSON
        with open(fname, 'w') as fh:
            json.dump(stmt_hashes, fh, indent=1)

        stmts_by_hash = {stmt.get_hash(): stmt for stmt in stmts}

        # Now do sanity checks
        preassembled_stmts_from_hashes = copy.deepcopy([all_stmts_by_hash[h]
                                                       for h in stmt_hashes])
        for stmt in preassembled_stmts_from_hashes:
            stmt.evidence = [e for e in stmt.evidence
                             if e.source_api == reader]

            assert len(stmt.evidence) == len(stmts_by_hash[stmt.get_hash()].evidence)
            assert set([e.source_hash for e in stmt.evidence]) == \
                set([e.source_hash for e in stmts_by_hash[stmt.get_hash()].evidence])