import numpy
import pickle
from collections import defaultdict, Counter
import scipy.optimize
import matplotlib.pyplot as plt
from texttable import Texttable
from indra_db.client.principal.curation import get_curations
from indra_db import get_primary_db


db = get_primary_db()


def belief(num_ev, pr, ps):
    return 1 - (ps + pr ** num_ev)


def logl(correct_by_num_ev, pr, ps):
    """Return log likelihood of belief model parameters given data."""
    ll = 0
    for num_ev, corrects in correct_by_num_ev.items():
        ll += sum((numpy.log(belief(num_ev, pr, ps)) if c else
                   numpy.log(1-belief(num_ev, pr, ps)))
                  for c in corrects)
    return ll


def optimize(correct_by_num_ev):
    """Return maximum likelihood parameters for belief model."""
    # The function being optimized is the negative log likelihood
    fun = lambda x: -logl(correct_by_num_ev, x[0], x[1])
    # Both parameters have to be between 0 and 1
    bounds = [(0, 1), (0, 1)]
    # 1 - ps - pr cannot be negative
    constr = lambda x: (1-x[0]-x[1])
    res = scipy.optimize.minimize(fun,
                                  x0=[0.05, 0.3],  # Initial guess: default
                                  bounds=bounds,
                                  constraints=[{'type': 'ineq',
                                                'fun': constr}])
    return res.x


def plot_curations(sources):
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
            curations[cur.pa_hash][cur.source_hash].append(cur.tag)

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
            print('Not enough curations for %s from %s: %s' %
                  (stmts_dict[pa_hash].uuid, source, stmts_dict[pa_hash]))
            continue
        # We can now assign 0 or 1 to each evidence's curation(s), resolve
        # any inconsistencies at the level of a single evidence.
        for source_hash, ev_curs in stmt_curs.items():
            corrects = [1 if cur in ('correct', 'hypothesis', 'act_vs_amt')
                        else 0 for cur in ev_curs]
            if any(corrects) and not all(corrects):
                print('Suspicious curation: (%s, %s), %s. Assuming overall'
                      ' incorrect.' % (pa_hash, source_hash, str(ev_curs)))
            overall_cur = 1 if all(corrects) else 0
            # We also need to make sure that if the same evidence hash appears
            # multiple times, we count it the right number of times
            overall_cur_by_num_ev_hash = \
                [overall_cur] * ev_hash_count[source_hash]
            full_curations[pa_hash] += overall_cur_by_num_ev_hash

    # Now we aggregate evidence-level correctness at the statement level and
    # assign 1 or 0 to the statement depending on whether any of its evidences
    # are correct or not. We also account for the case where the same Statement
    # was sampled multiple times.
    correct_by_num_ev = defaultdict(list)
    for pa_hash, corrects in full_curations.items():
        any_correct = 1 if any(corrects) else 0
        any_correct_by_num_sampled = [any_correct] * stmt_counts[pa_hash]
        correct_by_num_ev[len(corrects)] += any_correct_by_num_sampled

    opt_r, opt_s = optimize(correct_by_num_ev)
    print('Maximum likelihood random error: %.3f' % opt_r)
    print('Maximum likelihood systematic error: %.3f' % opt_s)

    # Finally, calculate the mean of correctness by number of evidence
    num_evs = sorted(correct_by_num_ev.keys())
    means = [numpy.mean(correct_by_num_ev[n]) for n in num_evs]
    # Stderr of proportion is sqrt(pq/n)
    std = [numpy.sqrt((numpy.mean(correct_by_num_ev[n]) *
                            (1 - numpy.mean(correct_by_num_ev[n]))) /
                            len(correct_by_num_ev[n]))
                for n in num_evs]
    beliefs = [belief(n, opt_r, opt_s) for n in num_evs]

    # Print table of results before plotting
    table = Texttable()
    table_data = [['Num Evs', 'Count', 'Num Correct', 'Pct', 'Std']]
    for i, num_ev in enumerate(num_evs):
        table_row = [num_ev, len(correct_by_num_ev[num_ev]),
                     sum(correct_by_num_ev[num_ev]), means[i], std[i]]
        table_data.append(table_row)
    table.add_rows(table_data)
    print(table.draw())

    plt.errorbar(num_evs, means, yerr=std, fmt='bo-',
                 label='Empirical mean correctness')
    plt.plot(num_evs, beliefs, 'ro-', label='Optimized belief')
    plt.ylim(0, 1)
    plt.grid(True)
    plt.xticks(num_evs)
    plt.xlabel('Number of evidence per INDRA Statement')
    plt.legend(loc='lower right')
    plt.show()


if __name__ == '__main__':
    input_source = 'reach'
    with open('../../data/curation/bioexp_%s_sample_uncurated.pkl'
              % input_source, 'rb') as fh:
        stmts = pickle.load(fh)
    with open('../../data/curation/bioexp_%s_sample_tsv.pkl' % input_source,
              'rb') as fh:
        tsv_stmts = pickle.load(fh)
        for stmt in tsv_stmts:
            stmt.evidence = [e for e in stmt.evidence
                             if e.source_api == input_source]
            stmts.append(stmt)

    stmts_dict = {stmt.get_hash(): stmt for stmt in stmts}
    stmt_counts = Counter(stmt.get_hash() for stmt in stmts)
    # Load all curations from the DB
    curation_sources = [('bioexp_paper_tsv', 'bioexp_paper_reach')]
    for source in curation_sources:
        plot_curations(source)
