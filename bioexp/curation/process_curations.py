import sys
import csv
import numpy
from collections import defaultdict
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
    curations = defaultdict(lambda: defaultdict(list))
    for source in sources:
        db_curations = get_curations(db, source=source)
        # For now, assume that all evidences for each statement have been
        # curated--otherwise we have to cross-reference back to the original
        # pickle to determine the number of evidences for the stmt/source
        # TODO Currently doesn't correctly handle cases of repeated sampling--
        # should probably use a nested dictionary
        for cur in db_curations:
            curations[cur.pa_hash][cur.source_hash].append(cur.tag)

    full_curations = defaultdict(list)
    for pa_hash, stmt_curs in curations.items():
        for source_hash, ev_curs in stmt_curs.items():
            corrects = [1 if cur in ('correct', 'hypothesis', 'act_vs_amt')
                        else 0 for cur in ev_curs]
            if any(corrects) and not all(corrects):
                print('Suspicious curation: (%s, %s), %s. Assuming overall'
                      ' incorrect.' % (pa_hash, source_hash, str(ev_curs)))
            overall_cur = 1 if all(corrects) else 0
            full_curations[pa_hash].append(overall_cur)
        # TODO: Cross-reference against assembly file to determine if all
        # curated
        # Filter to only curations where every entry for the
        # given UUID was curated
        #full_curations = {k: v for k, v in curations.items()
        #                  if all([vv is not None for vv in v])}
        # TODO: Cross-reference against sample file to determine if sampled
        # multiple times

    # Now organize the curations by number of evidence
    correct_by_num_ev = {}
    for pa_hash, corrects in full_curations.items():
        any_correct = 1 if any(corrects) else 0
        if len(corrects) in correct_by_num_ev:
            correct_by_num_ev[len(corrects)].append(any_correct)
        else:
            correct_by_num_ev[len(corrects)] = [any_correct]

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
    # Load all curations from the DB
    curation_sources = [('bioexp_paper_tsv', 'bioexp_paper_reach')]
    for source in curation_sources:
        plot_curations(source)
