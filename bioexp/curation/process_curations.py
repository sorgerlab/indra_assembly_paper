import json
import pickle
import logging
import itertools
from multiprocessing import Pool
from os.path import dirname, abspath, join
import corner
import scipy.optimize
from texttable import Texttable
import matplotlib.pyplot as plt
from collections import defaultdict, Counter
from indra_db import get_primary_db
from indra_db.client.principal.curation import get_curations
from bioexp.curation.belief_models import *
from bioexp.curation.model_fit import ModelFit, ens_sample


logger = logging.getLogger('process_curations')


db = get_primary_db()


def optimize(correct_by_num_ev):
    """Return maximum likelihood parameters for belief model."""
    # The function being optimized is the negative log likelihood
    fun = lambda x: -logl(correct_by_num_ev, x[0], x[1])
    # Both parameters have to be between 0 and 1
    bounds = [(0.01, 0.99), (0.01, 0.99)]
    res = scipy.optimize.minimize(fun,
                                  x0=[0.3, 0.05],  # Initial guess: default
                                  bounds=bounds)
    return res.x


def get_correctness_data(sources, stmts, aggregation='evidence'):
    stmts_dict = {stmt.get_hash(): stmt for stmt in stmts}
    stmt_counts = Counter(stmt.get_hash() for stmt in stmts)
    full_curations = get_full_curations(sources, stmts_dict,
                                        aggregation=aggregation)
    correct_by_num_ev = defaultdict(list)
    for pa_hash, corrects in full_curations.items():
        num_correct = sum(corrects)
        num_correct_by_num_sampled = [num_correct] * stmt_counts[pa_hash]
        correct_by_num_ev[len(corrects)] += num_correct_by_num_sampled
    return correct_by_num_ev


def get_raw_curations(sources, stmts_dict):
    # Curations are in a dict keyed by pa_hash and then by evidence source hash
    curations = defaultdict(lambda: defaultdict(list))
    # Iterate over all the curation sources
    for source in sources:
        # Get curations from DB from given curation source
        db_curations = get_curations(db, source=source)
        # We populate the curations dict with entries from the DB
        for cur in db_curations:
            if cur.pa_hash not in stmts_dict:
                print('Curation pa_hash is missing fron list of Statements '
                      'loaded from pickle: %s' %
                      str((cur.pa_hash, cur.source_hash, cur.tag, cur.source)))
            curations[cur.pa_hash][cur.source_hash].append(cur)
    return curations


def get_full_curations(sources, stmts_dict, aggregation='evidence'):
    curations = get_raw_curations(sources, stmts_dict)
    # Next we construct a dict of all curations that are "full" in that all
    # evidences of a given statement were curated, keyed by pa_hash
    full_curations = defaultdict(list)
    for pa_hash, stmt_curs in curations.items():
        # We need to make sure that all the evidence hashes were covered by the
        # curations in the DB. Note that we cannot go by number of curations
        # since two subtly different evidences can have the same hash, and
        # multiple curations sometimes exist for the same evidence.
        ev_hashes = [e.get_source_hash() for e in stmts_dict[pa_hash].evidence]
        ev_hash_count = Counter(ev_hashes)
        if set(stmt_curs.keys()) != set(ev_hashes):
            # If not all evidences are covered by curations, we print enough
            # details to identify the statement to complete its curations.
            print('Not enough curations for %s: %s' %
                  (stmts_dict[pa_hash].uuid, stmts_dict[pa_hash]))
            continue
        # We can now assign 0 or 1 to each evidence's curation(s), resolve
        # any inconsistencies at the level of a single evidence.
        pmid_curations = defaultdict(list)
        for source_hash, ev_curs in stmt_curs.items():
            ev = _find_evidence_by_hash(stmts_dict[pa_hash], source_hash)
            corrects = [1 if cur.tag in ('correct', 'hypothesis', 'act_vs_amt')
                        else 0 for cur in ev_curs]
            if any(corrects) and not all(corrects):
                print('Suspicious curation: (%s, %s), %s. Assuming overall'
                      ' incorrect.' % (pa_hash, source_hash,
                                       str([(c.tag, c.curator)
                                            for c in ev_curs])))
            overall_cur = 1 if all(corrects) else 0
            # We also need to make sure that if the same evidence hash appears
            # multiple times, we count it the right number of times
            overall_cur_by_num_ev_hash = \
                [overall_cur] * ev_hash_count[source_hash]
            pmid_curations[ev.pmid] += overall_cur_by_num_ev_hash
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
            full_curations[pa_hash] = \
                list(itertools.chain(*pmid_curations.values()))

    return full_curations


def _find_evidence_by_hash(stmt, source_hash):
    for ev in stmt.evidence:
        if ev.get_source_hash() == source_hash:
            return ev


def optimize_params(correct_by_num_ev):
    opt_r, opt_s = optimize(correct_by_num_ev)
    print('Maximum likelihood random error: %.3f' % opt_r)
    print('Maximum likelihood systematic error: %.3f' % opt_s)
    return opt_r, opt_s


def load_reach_curated_stmts():
    logger.info('Loading REACH statement pickles')
    with open('../../data/curation/bioexp_reach_sample_uncurated_19-12-14.pkl',
              'rb') as fh:
        stmts = pickle.load(fh)
    with open('../../data/curation/bioexp_reach_sample_uncurated_20-02-19.pkl',
              'rb') as fh:
        stmts += pickle.load(fh)
    with open('../../data/curation/bioexp_reach_sample_tsv.pkl',
              'rb') as fh:
        tsv_stmts = pickle.load(fh)
        for stmt in tsv_stmts:
            stmt.evidence = [e for e in stmt.evidence
                             if e.source_api == 'reach']
            stmts.append(stmt)
    return stmts


if __name__ == '__main__':
    plt.ion()

    stmts = load_reach_curated_stmts()
    source_list = ('bioexp_paper_tsv', 'bioexp_paper_reach')
    ev_correct_by_num_ev = get_correctness_data(source_list, stmts,
                                                aggregation='evidence')
    ev_correct_by_num_pmid = get_correctness_data(source_list, stmts,
                                                aggregation='pmid')

    # Load evidence frequency data
    ev_dist_path = join(dirname(abspath(__file__)),
                        'reach_stmt_evidence_distribution.json')
    with open(ev_dist_path, 'rt') as f:
        ev_dist = json.load(f)
        # Convert string keys to integer keys
        ev_dist = {int(k): v for k, v in ev_dist.items()}

    # Load PMID frequency data
    pmid_dist_path = join(dirname(abspath(__file__)),
                          'reach_stmt_pmid_distribution.json')
    with open(pmid_dist_path, 'rt') as f:
        pmid_dist = json.load(f)
        # Convert string keys to integer keys
        pmid_dist = {int(k): v for k, v in pmid_dist.items()}

    aggregations = {'pmid': (ev_correct_by_num_pmid, pmid_dist),
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
            filename = f'{model_name}_sampler.pkl'
            print(f'Saving to {filename}')
            with open(filename, 'wb') as f:
                sampler.pool = None
                pickle.dump((mf, sampler), f)
            results.append((model_name, mf, sampler))
            mf.plot_ev_fit(sampler, model_name)
            mf.plot_stmt_fit(sampler, model_name)

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

