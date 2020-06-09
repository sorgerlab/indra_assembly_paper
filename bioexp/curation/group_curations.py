import itertools
from indra.tools import assemble_corpus as ac
from bioexp.util import pkldump
from process_curations import get_correctness_data, load_curated_pkl_files, \
                              get_full_curations


reader_input = {'reach': {
                 'pkl_list': [
                    'bioexp_reach_sample_uncurated_19-12-14.pkl',
                    'bioexp_reach_sample_uncurated_20-02-19.pkl',
                    'bioexp_reach_sample_tsv.pkl'],
                 'source_list': ['bioexp_paper_reach', 'bioexp_paper_tsv']},
               'rlimsp': {
                 'pkl_list': ['bioexp_rlimsp_sample_uncurated.pkl'],
                 'source_list': ['bioexp_paper_rlimsp']},
               'trips': {
                 'pkl_list': ['bioexp_trips_sample_uncurated.pkl'],
                 'source_list': ['bioexp_paper_trips']}
               }


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
# at least 1 extraction is correct; for the ones that are incorrect, check
# if there is evidence from other readers; if so, add this to a pickle
# to be curated
all_stmts = ac.load_statements('../../data/bioexp_asmb_preassembled.pkl')
all_stmts_by_hash = {stmt.get_hash(): stmt for stmt in all_stmts}

# Lists or sets? Sets may eliminate necessary multiplicity?
corr_hashes = set()
incorr_hashes_1_src = set()
incorr_hashes_multi_src = set()
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
                incorr_hashes_multi_src.add(pa_hash)
            else:
                incorr_hashes_1_src.add(pa_hash)

source_list = ['bioexp_paper_reach', 'bioexp_paper_trips',
               'bioexp_paper_rlimsp', 'bioexp_paper_multi']
curations['multi'] = get_full_curations(source_list, all_stmts_by_hash,
                                        filter_hashes=incorr_stmts_multi_src)
# At this point, if a curation is all zeros, we can confirm that it is
# incorrect
incorr_hashes_multi_curated = set()
for pa_hash, corrects in curations['multi'].items():
    if sum(corrects) > 0:
        corr_hashes.add(pa_hash)
    else:
        incorr_hashes_multi_curated.add(pa_hash)

incorr_hashes = incorr_hashes_1_src | incorr_hashes_multi_curated
curated = corr_hashes | incorr_hashes
uncurated = all_hashes - curated

curated_stmts = [all_stmts_by_hash[h] for h in curated]

# Group curated stmts into bins
bins = [0., 0.6, 0.8, 0.9, 0.95, 0.99, 1.0]
stmts_by_belief = {}
for bin_ix in range(len(bins) - 1):
    lb = bins[bin_ix] # Lower bound
    ub = bins[bin_ix + 1] # Upper bound
    bin_stmts = [s for s in curated_stmts
                 if s.belief > lb and s.belief <= ub]
    # Get correctness stats for statements in each bin
    n_corr = 0
    n_incorr = 0
    for stmt in bin_stmts:
        pa_hash = stmt.get_hash()
        if pa_hash in corr_hashes:
            n_corr += 1
        elif pa_hash in incorr_hashes:
            n_incorr += 1
        else:
            assert False
    n_total = n_corr + n_incorr
    pct_corr = n_corr / n_total
    stmts_by_belief[bin_ix] = {'lb': lb, 'ub': ub, 'stmts': bin_stmts,
                               'n_corr': n_corr, 'n_incorr': n_incorr,
                               'n_total': n_total, 'pct_corr': pct_corr}
    print(f'{lb}-{ub}: {n_corr} / {n_total} = {pct_corr}')

pkldump(stmts_by_belief, 'multi_src_results')

# Load curations for the incorr_multi_src statements to determine if they
# have been fully curated
# - load curations from all sources, but filter to those in the specific
#   set

# Load all curations across all pickles
# Don't even have to stratify by source?

# Need to allow partially curated statements during load as long
# as one of them is correct


# In future, will also need to get curations for the multi-src pickle
# and use those to file into the correct or incorrect bins
# So, correct:
# - those statements from single readers that have been marked correct
# - those statements from multi-src pickle marked correct
# Incorrect:
# - curated stmts from a single reader with no other sources and no correct evs
# - curated stmts from multi-src pickle with no correct evs
#ac.dump_statements(incorr_stmts_multi_src, 'bioexp_incorr_multi_src.pkl')


#curated_stmts_by_hash = {stmt.get_hash(): stmt
#                         for stmt in curated_stmts}
#all_stmts_by_uuid = {stmt.uuid: stmt for stmt in all_stmts}

#print("Num curated stmts (hash):", len(all_hashes))
#print("Num curated correct:", len(corr_hashes))
#print("Num 1-src incorrect:", len(incorr_hashes_1_src))
#print("Num multi-src incorrect (req curation):", len(incorr_stmts_multi_src))


