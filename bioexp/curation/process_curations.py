import emcee
import pickle
import corner
import logging
import numpy as np
import scipy.optimize
from texttable import Texttable
import matplotlib.pyplot as plt
from collections import defaultdict, Counter
from indra_db import get_primary_db
from indra_db.client.principal.curation import get_curations


logger = logging.getLogger('process_curations')


db = get_primary_db()


def belief(num_ev, pr, ps):
    b = 1 - (ps + (1-ps) * (pr ** num_ev))
    return b


def logl(correct_by_num_ev, pr, ps):
    """Return log likelihood of belief model parameters given data."""
    ll = 0
    for num_ev, corrects in correct_by_num_ev.items():
        ll += sum((np.log(belief(num_ev, pr, ps)) if c else
                   np.log(1-belief(num_ev, pr, ps)))
                  for c in corrects)
    return ll


def log_prob(p, correct_by_num_ev):
    pr, ps = p
    # Uniform prior over [0, 1]
    if pr < 0 or ps < 0 or pr > 1 or ps > 1:
        return -np.inf
    else:
        return logl(correct_by_num_ev, pr, ps)


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


def preprocess_data(sources, stmts):
    stmts_dict = {stmt.get_hash(): stmt for stmt in stmts}
    stmt_counts = Counter(stmt.get_hash() for stmt in stmts)
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
            corrects = [1 if cur.tag in ('correct', 'hypothesis', 'act_vs_amt')
                        else 0 for cur in ev_curs]
            if any(corrects) and not all(corrects):
                print('Suspicious curation: (%s, %s), %s. Assuming overall'
                      ' incorrect.' % (pa_hash, source_hash,
                                       str([(c.tag, c.curator) for c in ev_curs])))
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
        npmid = len({ev.pmid for ev in stmts_dict[pa_hash].evidence})
        any_correct = 1 if any(corrects) else 0
        any_correct_by_num_sampled = [any_correct] * stmt_counts[pa_hash]
        correct_by_num_ev[len(corrects)] += any_correct_by_num_sampled
    return correct_by_num_ev


def optimize_params(correct_by_num_ev):
    opt_r, opt_s = optimize(correct_by_num_ev)
    print('Maximum likelihood random error: %.3f' % opt_r)
    print('Maximum likelihood systematic error: %.3f' % opt_s)
    return opt_r, opt_s


def plot_curations(correct_by_num_ev, opt_r, opt_s):
    # Finally, calculate the mean of correctness by number of evidence
    num_evs = sorted(correct_by_num_ev.keys())
    means = [np.mean(correct_by_num_ev[n]) for n in num_evs]
    # Stderr of proportion is sqrt(pq/n)
    std = [2*np.sqrt((np.mean(correct_by_num_ev[n]) *
                      (1 - np.mean(correct_by_num_ev[n]))) /
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


def load_reach_curated_stmts():
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
    return stmts


def get_posterior_samples(correct_by_num_ev, ndim=2, nwalkers=50,
                          nsteps=2000, nburn=100):
    logger.info('Sampling with %d walkers for %d steps' % (nwalkers, nsteps))
    # Initialize walkers across the interval [0, 1)
    p0 = np.random.rand(nwalkers, ndim)

    from multiprocessing import Pool
    with Pool() as pool:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob,
                                        args=[correct_by_num_ev],
                                        pool=pool)
        sampler.run_mcmc(p0, nsteps)
    logger.info('Finished sampling')
    return sampler.flatchain[nburn*nwalkers:]


def plot_posterior_samples(samples):
    # Plot the posterior parameter distribution
    plt.figure()
    corner.corner(samples, labels=['Rand.', 'Syst'])

    # Plot a few representative belief curves from the posterior
    num_evs = sorted(correct_by_num_ev.keys())
    for pr, ps in samples[:100]:
        beliefs = [belief(n, pr, ps) for n in num_evs]
        plt.plot(num_evs, beliefs, 'g-', alpha=0.1)


if __name__ == '__main__':
    plt.ion()
    stmts = load_reach_curated_stmts()
    source_list = ('bioexp_paper_tsv', 'bioexp_paper_reach')
    correct_by_num_ev = preprocess_data(source_list, stmts)

    # Basic optimization and max-likelihood estimates
    opt_r, opt_s = optimize_params(correct_by_num_ev)
    plot_curations(correct_by_num_ev, opt_r, opt_s)

    # Bayesian parameter inference
    samples = get_posterior_samples(correct_by_num_ev,
                                    ndim=2, nwalkers=10,
                                    nsteps=10000, nburn=100)

    plot_posterior_samples(samples)