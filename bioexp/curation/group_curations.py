import sys
import pickle
from os.path import abspath, dirname, join
from collections import Counter, defaultdict
from indra.tools import assemble_corpus as ac
from bioexp.curation.process_curations import \
            get_correctness_data, load_curated_pkl_files, get_full_curations, \
            reader_input, get_raw_curations, CURATIONS


def get_multi_reader_curations(reader_curations, reader_input,
                               all_stmts_by_hash):
    # For every curation hash across the pickles, set aside the ones where
    # at least 1 extraction is correct; for the ones that are incorrect, check
    # if there is evidence from other readers; if so, add this to a pickle
    # to be curated
    # TODO: Lists or sets? Sets may eliminate necessary multiplicity?
    corr_hashes = set()
    incorr_hashes_1_src = set()
    incorr_hashes_multi_src = set()
    all_hashes = set()
    for reader, rdr_curs in reader_curations.items():
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

    # Get curations from all DB curation tags for individual readers as
    # well as the multi-reader curations
    source_list = [source for reader_dict in reader_input.values()
                          for source in reader_dict['source_list']]
    source_list.append('bioexp_paper_multi')
    multi_curations = get_full_curations(source_list, all_stmts_by_hash,
                                         filter_hashes=incorr_hashes_multi_src)
    # At this point, if a curation is all zeros, we can confirm that it is
    # incorrect
    incorr_hashes_multi_curated = set()
    for pa_hash, corrects in multi_curations.items():
        if sum(corrects) > 0:
            corr_hashes.add(pa_hash)
        else:
            incorr_hashes_multi_curated.add(pa_hash)

    # Put together the results
    incorr_hashes = incorr_hashes_1_src | incorr_hashes_multi_curated
    curated = corr_hashes | incorr_hashes
    uncurated = all_hashes - curated
    results = {'multi_src_curations': multi_curations,
               'correct_hashes': corr_hashes,
               'incorrect_hashes': incorr_hashes,
               'multi_src': incorr_hashes_multi_src,
               'uncurated': uncurated}
    return results


def get_single_reader_curations(curations):
    # Version 2.0: Take all curated hashes (i.e., statements for which we can
    # establish correct or incorrect definitively) and strip out the other
    # evidences
    # curations is a set dict of PA hashes with the vector of 0/1 correctness
    # scores for evidences
    corr_hashes = set()
    incorr_hashes = set()
    for pa_hash, corrects in curations.items():
        if sum(corrects) > 0:
            corr_hashes.add(pa_hash)
        else:
            incorr_hashes.add(pa_hash)
    results = {'correct_hashes': corr_hashes,
               'incorrect_hashes': incorr_hashes}
    return results


# Load curations for all readers
def load_reader_curations(reader_input):
    """Get a dict of curations for each reader."""
    curations = {}
    stmts_by_reader = {}
    for reader, rd_dict in reader_input.items():
        reader_stmts = load_curated_pkl_files(rd_dict['pkl_list'])
        reader_stmts_dict = {stmt.get_hash(): stmt for stmt in reader_stmts}
        #curations[reader] = get_correctness_data(rd_dict['source_list'],
        #                           reader_stmts, aggregation='evidence')
        curations[reader] = get_full_curations(rd_dict['source_list'],
                                 reader_stmts_dict, aggregation='evidence',
                                 allow_incomplete_correct=True)
    return curations


def dump_dataset(stmts_by_hash, corr_hashes, incorr_hashes, filename):
    # Prepare dataset for statistical modeling
    kge_data = []
    # For each stmt get the evidence from each source
    curated_hashes = corr_hashes | incorr_hashes
    for ix, stmt_hash in enumerate(curated_hashes):
        # Get the statement
        stmt = stmts_by_hash[stmt_hash]
        # Get the number of evidences for each source
        sources = [ev.source_api for ev in stmt.evidence]
        source_entry = dict(Counter(sources))
        # Get the overall correctness status from the multi_results dict
        corr = 1 if stmt.get_hash() in corr_hashes else 0
        source_entry['correct'] = corr
        source_entry['stmt_type'] = stmt.__class__.__name__
        # Add basic statement data (useful for linking to knowledge graph
        # embedding-based link predict)
        agent_names = [ag.name for ag in stmt.agent_list() if ag is not None]
        # Complex > 3, translocations, autophosphorylations will be skipped
        if len(agent_names) == 2:
            kge_entry = {'stmt_num': ix,
                         'stmt_hash': stmt.get_hash(),
                         'agA_name': agent_names[0],
                         'stmt_type': stmt.__class__.__name__,
                         'agB_name': agent_names[1],
                         'correct': corr}
            kge_entry.update(source_entry)
            kge_data.append(kge_entry)

    with open(filename, 'wb') as f:
        pickle.dump(kge_data, f)

    return kge_data


def get_combined_curations(source_list, stmts_by_hash, filename,
                           add_supports=False, allow_incomplete=False,
                           allow_incomplete_correct=True):
    """This function creates a custom curation data structure to
     facilitate downstream analysis of the curation data. The returned
     data structure is a list with dict entries. Each dict corresponds
     to a curated statement and carries necessary metadata about the
     statement itself, sources supporting the statement, as well as the
     overall correcness of the statement per all aggragated curations."""
    # Prepare dataset for statistical modeling
    cur_data = []
    # Get curations for all sources
    full_curations = get_full_curations(source_list, stmts_by_hash,
                                        allow_incomplete=allow_incomplete,
                                        allow_incomplete_correct=allow_incomplete_correct)
    # Build up a dictionary of hashes associated with each curation tag
    stmt_tags_by_hash = defaultdict(list)
    for source_tag in source_list:
        db_curations = get_raw_curations([source_tag], stmts_by_hash)
        db_cur_hashes = list(db_curations.keys())
        for cur_hash in db_cur_hashes:
            stmt_tags_by_hash[cur_hash].append(source_tag)
    stmt_tags_by_hash = dict(stmt_tags_by_hash)

    for ix, (stmt_hash, corrects) in enumerate(full_curations.items()):
        # Get the statement
        stmt = stmts_by_hash[stmt_hash]
        # Complex > 3, translocations, autophosphorylations will be skipped
        agents = stmt.real_agent_list()
        if len(agents) != 2:
            continue
        # Get the number of evidences for each source
        sources = [ev.source_api for ev in stmt.evidence]
        source_entry = dict(Counter(sources))
        # Get the overall correctness status
        num_correct = sum(corrects)
        corr = 1 if num_correct > 0 else 0
        source_entry['correct'] = corr
        # Add in source evidence counts from supports
        if add_supports:
            source_entry['num_supports'] = len(stmt.supports)
            for supp_stmt in stmt.supports:
                supp_sources = [ev.source_api for ev in supp_stmt.evidence]
                supp_source_ctr = dict(Counter(supp_sources))
                # Add the supporting stmt sources to the stmt's own sources
                for source_api, source_ct in supp_source_ctr.items():
                    if source_api not in source_entry:
                        source_entry[source_api] = 0
                    source_entry[source_api] += source_ct
        # Add basic statement data (useful for linking to knowledge graph
        # embedding-based link predict)
        agA_ns, agA_id = agents[0].get_grounding()
        agA_name = agents[0].name
        agB_ns, agB_id = stmt.agent_list()[1].get_grounding()
        agB_name = agents[1].name

        cur_entry = {'stmt_num': ix,
                     'stmt_hash': stmt.get_hash(),
                     'agA_name': agA_name,
                     'agA_ns': agA_ns,
                     'agA_id': agA_id,
                     'stmt_type': stmt.__class__.__name__,
                     'agB_name': agB_name,
                     'agB_ns': agB_ns,
                     'agB_id': agB_id,
                     'correct': corr}
        # Add columns indicating which curation samples included this stmt
        for source_tag in source_list:
            has_tag = int(source_tag in stmt_tags_by_hash[stmt_hash])
            cur_entry[source_tag] = has_tag

        # Add this entry to the dataset
        cur_entry.update(source_entry)
        cur_data.append(cur_entry)

    with open(filename, 'wb') as f:
        pickle.dump(cur_data, f)

    return cur_data


def print_curation_stats(fname):
    with open(fname, 'rb') as fh:
        curs = pickle.load(fh)
    stmt_hashes = {c['stmt_hash'] for c in curs}
    mention_hashes = {c['source_hash'] for c in CURATIONS
                      if c['pa_hash'] in stmt_hashes}
    print('For %s, the number of unique Statements is %d and the number '
          'of corresponding mentions with curations is %d.' %
          (fname, len(stmt_hashes), len(mention_hashes)))


if __name__ == '__main__':
    # Prevent issues in pickling the results
    sys.setrecursionlimit(50000)

    output_dir = sys.argv[1]

    # Load the pickle file with all assembled statements
    asmb_pkl = join(dirname(abspath(__file__)), '..', '..', 'data',
                    'bioexp_asmb_preassembled.pkl')
    all_stmts = ac.load_statements(asmb_pkl)
    all_stmts_by_hash = {stmt.get_hash(): stmt for stmt in all_stmts}

    all_sources = [source for rdr, rdr_dict in reader_input.items()
                   for source in rdr_dict['source_list']]
    all_sources.append('bioexp_paper_multi')

    multireader_curation_dataset = get_combined_curations(
          all_sources, all_stmts_by_hash,
          join(output_dir, 'multireader_curation_dataset.pkl'),
          allow_incomplete=True)
    extended_curation_dataset = get_combined_curations(
          all_sources + ['bioexp_biogrid', 'bioexp_psp'],
          all_stmts_by_hash,
          join(output_dir, 'extended_curation_dataset.pkl'),
          allow_incomplete=True)
