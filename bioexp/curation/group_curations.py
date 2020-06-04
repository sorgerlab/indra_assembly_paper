import itertools
from indra.tools import assemble_corpus as ac
from process_curations import get_correctness_data, load_curated_pkl_files, \
                              get_full_curations


reader_input = {'reach': {
                 'pkl_list': [
                    'bioexp_reach_sample_uncurated_19-12-14.pkl',
                    'bioexp_reach_sample_uncurated_20-02-19.pkl',
                    'bioexp_reach_sample_tsv.pkl'],
                 'source_list': ['bioexp_paper_reach',
                                 'bioexp_paper_tsv']},
               'rlimsp': {
                 'pkl_list': ['bioexp_rlimsp_sample_uncurated.pkl'],
                 'source_list': ['bioexp_paper_rlimsp']},
               'trips': {
                 'pkl_list': ['bioexp_trips_sample_uncurated.pkl'],
                 'source_list': ['bioexp_paper_trips']}
               }

# Skipped the following:
# 'bioexp_paper_tsv',

curated_stmts = []
curations = {}
for reader, rd_dict in reader_input.items():
    stmts = load_curated_pkl_files(rd_dict['pkl_list'])
    curated_stmts.extend(stmts)
    stmts_dict = {stmt.get_hash(): stmt for stmt in stmts}
    curations[reader] = get_full_curations(rd_dict['source_list'], stmts_dict,
                                   aggregation='evidence')

for r1, r2 in itertools.combinations(curations.keys(), 2):
    r1_cur = set(curations[r1].keys())
    r2_cur = set(curations[r2].keys())
    overlap_hashes = r1_cur.intersection(r2_cur)
    print(r1, r2, ": overlap ", len(overlap_hashes), " hashes")

# For every curation hash across the 3 pickles, set aside the ones where
# at least extraction is correct; for the ones that are incorrect, check
# if there is evidence from other readers; if so, add this to a pickle
# to be curated
all_stmts = ac.load_statements('../../data/bioexp_asmb_preassembled.pkl')
all_stmts_by_hash = {stmt.get_hash(): stmt for stmt in all_stmts}
corr_hashes = set()
incorr_hashes_1_src = set()
incorr_stmts_multi_src = []
not_found = set()
all_hashes = set()
for reader, rdr_curs in curations.items():
    for pa_hash, corrects in rdr_curs.items():
        all_hashes.add(pa_hash)
        if sum(corrects) > 0:
            corr_hashes.add(pa_hash)
        else:
            # Get the statement from the assembled pickle
            stmt = all_stmts_by_hash[pa_hash]
            sources = set([ev.source_api for ev in stmt.evidence])
            if len(sources) > 1:
                incorr_stmts_multi_src.append(stmt)
            else:
                incorr_hashes_1_src.add(pa_hash)

# In future, will also need to get curations for the multi-src pickle
# and use those to file into the correct or incorrect bins
# So, correct:
# - those statements from single readers that have been marked correct
# - those statements from multi-src pickle marked correct
# Incorrect:
# - curated stmts from a single reader with no other sources and no correct evs
# - curated stmts from multi-src pickle with no correct evs
ac.dump_statements(incorr_stmts_multi_src, 'bioexp_incorr_multi_src.pkl')


curated_stmts_by_hash = {stmt.get_hash(): stmt
                         for stmt in curated_stmts}
all_stmts_by_uuid = {stmt.uuid: stmt for stmt in all_stmts}

print("Num curated stmts (hash):", len(all_hashes))
print("Num curated correct:", len(corr_hashes))
print("Num 1-src incorrect:", len(incorr_hashes_1_src))
print("Num multi-src incorrect (req curation):", len(incorr_stmts_multi_src))


