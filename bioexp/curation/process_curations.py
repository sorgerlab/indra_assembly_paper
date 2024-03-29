import copy
import sys
import json
import pickle
import logging
import itertools
from multiprocessing import Pool
from os import pardir
from os.path import dirname, abspath, join
from texttable import Texttable
import matplotlib.pyplot as plt
from collections import defaultdict, Counter
from indra.tools import assemble_corpus as ac
from bioexp.util import prefixed_file, pkldump
from bioexp.curation.belief_models import *
from bioexp.curation.model_fit import ModelFit, ens_sample

logger = logging.getLogger('process_curations')
here = dirname(abspath(__file__))
curation_data = join(here, pardir, pardir, 'data', 'curation')


def dist_path(reader, dist_type):
    return prefixed_file(f'{reader}_stmt_{dist_type}_distribution', 'json')


def read_curations():
    cur_fname = join(curation_data, 'indra_assembly_curations.json')
    with open(cur_fname, 'r') as fh:
        return json.load(fh)


CURATIONS = read_curations()


reader_input = {
   'reach': {
     'pkl_list': [
        'bioexp_reach_sample_uncurated_19-12-14.pkl',
        'bioexp_reach_sample_uncurated_20-02-19.pkl',
        'bioexp_reach_sample_tsv.pkl'],
     'source_list': ['bioexp_paper_reach', 'bioexp_paper_tsv'],
     'belief_model': 'reach_orig_belief_stmt_evidence_sampler',
     'ev_dist_path': dist_path('reach', 'evidence'),
     'pmid_dist_path': dist_path('reach', 'pmid'),
    },
   'rlimsp': {
     'pkl_list': ['bioexp_rlimsp_sample_uncurated.pkl'],
     'source_list': ['bioexp_paper_rlimsp'],
     'belief_model': 'rlimsp_orig_belief_stmt_evidence_sampler',
     'ev_dist_path': dist_path('rlimsp', 'evidence'),
     'pmid_dist_path': dist_path('rlimsp', 'pmid'),
    },
   'trips': {
     'pkl_list': ['bioexp_trips_sample_uncurated.pkl'],
     'source_list': ['bioexp_paper_trips'],
     'belief_model': 'trips_orig_belief_stmt_evidence_sampler',
     'ev_dist_path': dist_path('trips', 'evidence'),
     'pmid_dist_path': dist_path('trips', 'pmid'),
    },
   'sparser': {
     'pkl_list': ['bioexp_sparser_sample_uncurated.pkl'],
     'source_list': ['bioexp_paper_sparser', 'bioexp_sparser'],
     'belief_model': 'sparser_orig_belief_stmt_evidence_sampler',
     'ev_dist_path': dist_path('sparser', 'evidence'),
     'pmid_dist_path': dist_path('sparser', 'pmid'),
    },
   'medscan': {
     'pkl_list': ['bioexp_medscan_sample_uncurated.pkl'],
     'source_list': ['bioexp_paper_medscan', 'bioexp_medscan'],
     'belief_model': 'medscan_orig_belief_stmt_evidence_sampler',
     'ev_dist_path': dist_path('medscan', 'evidence'),
     'pmid_dist_path': dist_path('medscan', 'pmid'),
    },
}


def get_curations_for_reader(reader, all_stmts, aggregation, **kwargs):
    """Get correctness data for a given reader based on reader_input info.

    Parameters
    ----------
    reader : str
        Name of the reader, e.g. "reach".
    all_stmts : list[indra.statements.Statement]
        A list of all statements in the assembled corpus.
    aggregation: str
        'evidence' to aggregate by distinct evidences, 'pmid' to
        aggregate by distinct PMIDs.

    Returns
    -------
    Dict
        Dictionary mapping total evidence (or PMID) counts to a list of the
        number of correct evidences/PMIDs for statements with that total
        evidence/PMID count.
    """
    # Get the info for this reader
    if reader in reader_input:
        ri = reader_input[reader]
        pkl_list = ri['pkl_list']
        source_list = ri['source_list']
    else:
        raise ValueError("Reader %s not supported." % reader)

    stmts = load_curated_pkl_files(pkl_list, all_stmts, reader)
    ev_corrects = get_correctness_data(source_list, stmts,
                                       aggregation=aggregation, **kwargs)
    return ev_corrects

# For debugging purposes, note module level global variables
ALL_HASHES = set()
ALL_MENTIONS = 0


def get_correctness_data(sources, stmts, aggregation='evidence',
                         allow_incomplete=False,
                         allow_incomplete_correct=False):
    global ALL_HASHES
    global ALL_MENTIONS
    stmts_dict = {stmt.get_hash(): stmt for stmt in stmts}
    stmt_counts = Counter(stmt.get_hash() for stmt in stmts)
    full_curations = get_full_curations(sources, stmts_dict,
                            aggregation=aggregation,
                            allow_incomplete=allow_incomplete,
                            allow_incomplete_correct=allow_incomplete_correct)
    correct_by_num_ev = {}
    ALL_HASHES |= set(full_curations)
    for pa_hash, corrects in full_curations.items():
        stmt = stmts_dict[pa_hash]
        ALL_MENTIONS += len(stmt.evidence)
        num_correct = sum(corrects)
        num_correct_by_num_sampled = [num_correct] * stmt_counts[pa_hash]
        if len(stmt.evidence) not in correct_by_num_ev:
            correct_by_num_ev[len(stmt.evidence)] = num_correct_by_num_sampled
        else:
            correct_by_num_ev[len(stmt.evidence)] += num_correct_by_num_sampled
    # Print ALL_HASHES and ALL_MENTIONS here to get the total numbers
    logger.info("Total hashes up to this point: %d" % len(ALL_HASHES))
    logger.info("Total mentions up to this point: %d" % ALL_MENTIONS)
    return correct_by_num_ev


def get_full_curations(sources, stmts_dict, aggregation='evidence',
                       filter_hashes=None,
                       allow_incomplete=False,
                       allow_incomplete_correct=False):
    """This function converts raw curations organized by statement/evidence
    and applies a set of policies to determine correctness for each evidence.
    It then returns a dictionary mapping statement hashes to lists of
    correctness values (0 for incorrect, 1 for correct) for each evidence in
    the statement."""

    curations = get_raw_curations(sources, stmts_dict)
    # Next we construct a dict of all curations that are "full" in that all
    # evidences of a given statement were curated, keyed by pa_hash
    full_curations = defaultdict(list)
    for pa_hash, stmt_curs in curations.items():
        if filter_hashes and pa_hash not in filter_hashes:
            continue
        if pa_hash not in stmts_dict:
            continue
        cur_stmt = stmts_dict[pa_hash]
        # We need to make sure that all the evidence hashes were covered by the
        # curations. Note that we cannot go by number of curations
        # since two subtly different evidences can have the same hash, and
        # multiple curations sometimes exist for the same evidence.
        ev_hashes = [e.get_source_hash() for e in cur_stmt.evidence]
        ev_hash_count = Counter(ev_hashes)
        # We can now assign 0 or 1 to each evidence's curation(s), resolve
        # any inconsistencies at the level of a single evidence.
        pmid_curations = defaultdict(list)
        for source_hash, ev_curs in stmt_curs.items():
            ev = _find_evidence_by_hash(cur_stmt, source_hash)
            if not ev:
                continue
            corrects = [1 if cur['tag'] in
                                    ('correct', 'hypothesis', 'act_vs_amt')
                        else 0 for cur in ev_curs]
            if any(corrects) and not all(corrects):
                print('Suspicious curation: (%s, %s, %s), %s. Assuming overall'
                      ' incorrect.' % (pa_hash, source_hash, cur_stmt,
                                       str([(c['tag'], c['curator'])
                                            for c in ev_curs])))
            overall_cur = 1 if all(corrects) else 0
            # We also need to make sure that if the same evidence hash appears
            # multiple times, we count it the right number of times
            overall_cur_by_num_ev_hash = \
                [overall_cur] * ev_hash_count[source_hash]
            pmid_curations[ev.pmid] += overall_cur_by_num_ev_hash
        evidence_corrects = list(itertools.chain(*pmid_curations.values()))

        # The statement is complete, or if we don't care if it's not complete
        # let the statement through
        if allow_incomplete or set(stmt_curs.keys()) == set(ev_hashes):
            pass
        # The statement is incomplete but correct and we're allowing
        # incomplete correct stmts
        elif allow_incomplete_correct and any(evidence_corrects) and \
                    set(stmt_curs.keys()) != set(ev_hashes):
            print("Allowing incompletely curated correct stmt",
                  cur_stmt.uuid, cur_stmt)
        # Otherwise, the statement is incomplete AND either (so far) incorrect
        # or we're not allowing incomplete corrects
        else:
            # If not all evidences are covered by curations, we print enough
            # details to identify the statement to complete its curations.
            print('Not enough curations for %s: %s' %
                  (cur_stmt.uuid, cur_stmt))
            continue

        # If we aggregate at the pmid level, we determine correctness at the
        # level of each PMID to create the list of curations (whose number
        # will be equal to the number of PMIDs)
        if aggregation == 'pmid':
            full_curations[pa_hash] = \
                [1 if any(pmid_curs) else 0 for pmid_curs
                 in pmid_curations.values()]
        # Otherwise we simply take the list of evidence-level curations
        # i.e., we flatten the PMID-based curations into a single flat list.
        else:
            full_curations[pa_hash] = evidence_corrects
    return full_curations


def get_raw_curations(sources, stmts_dict):
    """Get curations associated with the given source tags.

    For each statement hash and evidence hash, this returns the list of
    corresponding curation dicts as is without any aggregation.

    Parameters
    ----------
    sources : list of str
        List of curation source tags (e.g. names of statement samples) to
        query the DB curations table with.
    stmts_dict : list of INDRA Statements
        Statements corresponding to a superset of the curated statements with
        the given tags.

    Returns
    -------
    dict
        Nested dictionary with the pa_hash as the first-level key, the
        source_hash as the second-level key, and the statement curation
        data itself as the innermost dictionary.
    """
    # Curations are in a dict keyed by pa_hash and then by evidence source hash
    curations = defaultdict(lambda: defaultdict(list))
    # Iterate over all the curation sources
    for source in sources:
        # Get curations from DB from given curation source
        # We populate the curations dict with entries from the DB
        db_curations = [cur for cur in CURATIONS if
                        cur['source'] == source]
        for cur in db_curations:
            if cur['pa_hash'] not in stmts_dict:
                print('Curation pa_hash is missing from list of Statements '
                      'loaded from pickle: %s' %
                      str((cur['pa_hash'], cur['source_hash'],
                           cur['tag'], cur['source'])))
            curations[cur['pa_hash']][cur['source_hash']].append(cur)
    logger.info('Loaded %d raw curations for sources %s' %
                (len(curations), str(sources)))
    return curations


def _find_evidence_by_hash(stmt, source_hash):
    for ev in stmt.evidence:
        if ev.get_source_hash() == source_hash:
            return ev


def load_curated_pkl_files(pkl_list, all_stmts, reader, use_jsons=True):
    """Load statements sampled for curation from a list of pickle files.

    Note that since the sampled pickles are large and not in version control,
    instead we load the statement hashes from the JSONs that are in version
    control and filter the preassembled statements for the hashes to create
    the same set of statements as the ones in the pickle files, unless
    use_jsons is set to False.
    """
    all_stmts_by_hash = {stmt.get_hash(): stmt for stmt in all_stmts}
    logger.info('Loading curation statement pickles')
    stmts = []
    for pkl_file in pkl_list:
        pkl_path = join(curation_data, pkl_file)
        if not use_jsons:
            logger.info('Loading %s' % pkl_path)
            with open(pkl_path, 'rb') as fh:
                pkl_stmts = pickle.load(fh)
                # Special handling for the pickle file for the TSV
                if pkl_file == 'bioexp_reach_sample_tsv.pkl':
                    for stmt in pkl_stmts:
                        stmt.evidence = [e for e in stmt.evidence
                                         if e.source_api == 'reach']
        else:
            json_path = pkl_path.replace('.pkl', '_hashes.json')
            logger.info('Loading %s' % json_path)
            with open(json_path, 'r') as fh:
                hashes = json.load(fh)
            pkl_stmts = copy.deepcopy([all_stmts_by_hash[hash]
                                       for hash in hashes])
            for stmt in pkl_stmts:
                stmt.evidence = [e for e in stmt.evidence
                                 if e.source_api == reader]
        stmts.extend(pkl_stmts)
    return stmts


def load_stmt_evidence_distribution(reader):
    """Return a dict of empirical evidence count distributions in the corpus."""
    ev_file = prefixed_file(f'{reader}_stmt_evidence_distribution', 'json')
    with open(ev_file, 'r') as fh:
        ev_probs = json.load(fh)
        ev_probs = {int(k): v for k, v in ev_probs.items()}
    return ev_probs


# MAIN -----------------------------------------------------------------


if __name__ == '__main__':
    plt.ion()

    reader = sys.argv[1]
    output_dir = sys.argv[2]

    # Load the pickle file with all assembled statements
    asmb_pkl = join(dirname(abspath(__file__)), '..', '..', 'data',
                    'bioexp_asmb_preassembled.pkl')
    all_stmts = ac.load_statements(asmb_pkl)

    ev_correct_by_num_ev = get_curations_for_reader(
                                    reader, all_stmts, aggregation='evidence',
                                    allow_incomplete=False)
    #ev_correct_by_num_pmid = get_curations_for_reader(
    #                                reader, all_stmts, aggregation='pmid',
    #                                allow_incomplete=False)

    # -- Everything below is for model fitting! --
    # Load evidence frequency data
    with open(reader_input[reader]['ev_dist_path'], 'rt') as f:
        ev_dist = json.load(f)
        # Convert string keys to integer keys
        ev_dist = {int(k): v for k, v in ev_dist.items()}

    # Load PMID frequency data
    with open(reader_input[reader]['pmid_dist_path'], 'rt') as f:
        pmid_dist = json.load(f)
        # Convert string keys to integer keys
        pmid_dist = {int(k): v for k, v in pmid_dist.items()}

    aggregations = {#'pmid': (ev_correct_by_num_pmid, pmid_dist),}
                    'evidence': (ev_correct_by_num_ev, ev_dist)}
    models = {
        'orig_belief_ev': OrigBeliefEv,
        'orig_belief_stmt': OrigBeliefStmt,
        'binom_ev': BinomialEv,
        'binom_stmt': BinomialStmt,
        'betabinom_ev': BetaBinomialEv,
        'betabinom_stmt': BetaBinomialStmt
        }
    results = []
    for aggregation_type, (data, ev_dist_weights) in aggregations.items():
        for model_name, model_class in models.items():
            model_name = f'{model_name}_{aggregation_type}'
            model = model_class(weights=ev_dist_weights)
            print(f'Fitting {model_name}')
            mf = ModelFit(model, data)
            nwalkers, burn_steps, sample_steps = (100, 100, 100)
            with Pool() as pool:
                sampler = ens_sample(mf, nwalkers, burn_steps, sample_steps,
                                     pool=pool)
            filename = f'{reader}_{model_name}_sampler'
            sampler.pool = None
            pkldump((mf, sampler), filename)
            results.append((model_name, mf, sampler))
            mf.plot_ev_fit(sampler, model_name)
            mf.plot_stmt_fit(sampler, model_name, 'red')

    stmt_lkls = []
    stmt_lkls_wt = []
    labels = []
    for model_name, mf, sampler in results:
        labels.append(model_name)
        stmt_lkls.append(mf.stmt_err(sampler))
        stmt_lkls_wt.append(mf.stmt_err(sampler, ev_dist))
    plt.figure()
    plt.bar(range(len(stmt_lkls)), stmt_lkls, tick_label=labels)
    plt.bar(range(len(stmt_lkls_wt)), stmt_lkls_wt,
            tick_label=[f'{l}_wts' for l in labels])
    plt.ylim(bottom=250)

    # Print table of results
    table = Texttable()
    table_data = [('Model', '-log(Max Lkl)', '-log(Max Lkl) Wtd.')]
    table_data.extend(zip(labels, stmt_lkls, stmt_lkls_wt))
    table.add_rows(table_data)
    print(table.draw())

    # Pickle results
    results_path = join(output_dir, f'fig4_model_fit_results_{reader}.pkl')
    with open(results_path, 'wb') as f:
        pickle.dump(results, f)

