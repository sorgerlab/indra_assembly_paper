import emcee
import pickle
import corner
import logging
from multiprocessing import Pool
import numpy as np
import scipy.optimize
from scipy.stats import binom
from scipy.special import betaln
from texttable import Texttable
import matplotlib.pyplot as plt
from collections import defaultdict, Counter
from indra_db import get_primary_db
from indra_db.client.principal.curation import get_curations
from bioexp.curation.belief_models import *
from bioexp.curation.model_fit import ModelFit, ens_sample


logger = logging.getLogger('process_curations')


db = get_primary_db()



'''
def logl(correct_by_num_ev, pr, ps):
    """Return log likelihood of belief model parameters given data."""
    ll = 0
    for num_ev, corrects in correct_by_num_ev.items():
        ll += sum((np.log(belief(num_ev, pr, ps)) if c else
                   np.log(1-belief(num_ev, pr, ps)))
                  for c in corrects)
    return ll

def logl(correct_by_num_ev, pr, ps):
    """Return log likelihood of belief model parameters given data."""
    ll = 0
    for num_ev, num_corrects in correct_by_num_ev.items():
        for num_correct in num_corrects:
            if num_correct == 0:
                ll += np.log(ps + (1-ps)*binom.pmf(0, n=num_ev, p=1-pr))
            else:
                ll += np.log((1-ps) * binom.pmf(num_correct, n=num_ev, p=1-pr))
    return ll
'''


def logl(correct_by_num_ev, alpha, beta):
    """Return log likelihood of belief model parameters given data."""
    ll = 0
    for num_ev, num_corrects in correct_by_num_ev.items():
        n = num_ev
        for num_correct in num_corrects:
            k = num_correct
            b1 = betaln(k + alpha, n - k + beta)
            b2 = betaln(alpha, beta)
            nck = -betaln(1 + n - k, 1 + k) - np.log(n + 1)
            ll += nck + b1 - b2
    return ll




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


def get_statement_correctness_data(sources, stmts):
    stmts_dict = {stmt.get_hash(): stmt for stmt in stmts}
    stmt_counts = Counter(stmt.get_hash() for stmt in stmts)
    full_curations = get_full_curations(sources, stmts_dict)

    # Now we aggregate evidence-level correctness at the statement level and
    # assign 1 or 0 to the statement depending on whether any of its evidences
    # are correct or not. We also account for the case where the same Statement
    # was sampled multiple times.
    correct_by_num_ev = defaultdict(list)
    for pa_hash, corrects in full_curations.items():
        any_correct = 1 if any(corrects) else 0
        any_correct_by_num_sampled = [any_correct] * stmt_counts[pa_hash]
        correct_by_num_ev[len(corrects)] += any_correct_by_num_sampled
    return correct_by_num_ev


def get_evidence_correctness_data(sources, stmts):
    stmts_dict = {stmt.get_hash(): stmt for stmt in stmts}
    stmt_counts = Counter(stmt.get_hash() for stmt in stmts)
    full_curations = get_full_curations(sources, stmts_dict)
    correct_by_num_ev = defaultdict(list)
    for pa_hash, corrects in full_curations.items():
        num_correct = sum(corrects)
        num_correct_by_num_sampled = [num_correct] * stmt_counts[pa_hash]
        correct_by_num_ev[len(corrects)] += num_correct_by_num_sampled
    return correct_by_num_ev


def get_full_curations(sources, stmts_dict):
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
    return full_curations


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

    plt.errorbar(num_evs, means, yerr=std, fmt='bo-', ls='none',
                 label='Empirical mean correctness')
    plt.plot(num_evs, beliefs, 'ro-', label='Optimized belief')
    plt.ylim(0, 1)
    plt.grid(True)
    plt.xticks(num_evs)
    plt.xlabel('Number of evidence per INDRA Statement')
    plt.legend(loc='lower right')
    plt.show()


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


"""
def get_posterior_samples(correct_by_num_ev, ndim=2, nwalkers=50,
                          nsteps=2000, nburn=100, p0=None):
    logger.info('Sampling with %d walkers for %d steps' % (nwalkers, nsteps))
    # Initialize walkers across the interval [0, 1)
    if p0 is None:
        p0 = np.random.rand(nwalkers, ndim)

    from multiprocessing import Pool
    with Pool() as pool:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob,
                                        args=[correct_by_num_ev],
                                        pool=pool)
        sampler.run_mcmc(p0, nsteps)
    logger.info('Finished sampling')
    return sampler.flatchain[nburn*nwalkers:]
"""

def plot_posterior_samples(samples):
    # Plot the posterior parameter distribution
    plt.figure()
    corner.corner(samples, labels=['Rand.', 'Syst'])

    # Plot a few representative belief curves from the posterior
    num_evs = sorted(ev_correct_by_num_ev.keys())
    for pr, ps in samples[:100]:
        beliefs = [belief(n, pr, ps) for n in num_evs]
        plt.plot(num_evs, beliefs, 'g-', alpha=0.1)


if __name__ == '__main__':
    #with open('orig_belief_stmt_sampler.pkl', 'rb') as f:
    #    (mf, sampler) = pickle.load(f)

    plt.ion()
    """
    stmts = load_reach_curated_stmts()
    source_list = ('bioexp_paper_tsv', 'bioexp_paper_reach')
    stmt_correct_by_num_ev = get_statement_correctness_data(source_list, stmts)
    ev_correct_by_num_ev = get_evidence_correctness_data(source_list, stmts)
    """
    with open('correct_by_num_ev.pkl', 'rb') as f:
        ev_correct_by_num_ev = pickle.load(f)

    # Basic optimization and max-likelihood estimates
    #opt_r, opt_s = optimize_params(ev_correct_by_num_ev)
    #plot_curations(ev_correct_by_num_ev, opt_r, opt_s)

    # Bayesian parameter inference
    #samples = get_posterior_samples(ev_correct_by_num_ev,
    #                                ndim=2, nwalkers=10,
    #                                nsteps=10000, nburn=100)
    #plot_posterior_samples(samples)

    be = BinomialEv()
    bs = BinomialStmt()
    bbe = BetaBinomialEv()
    bbs = BetaBinomialStmt()
    obe = OrigBeliefEv()
    obs = OrigBeliefStmt()
    models = [('orig_belief_ev', obe), ('orig_belief_stmt', obs),
              ('binom_ev', be), ('binom_stmt', bs),
              ('betabinom_ev', be), ('betabinom_stmt', bbs)]
    #models = [('orig_belief_ev', obe), ('orig_belief_stmt', obs)]
              #('betabinom_ev', bbe), ('betabinom_stmt', bbs)]

    results = []
    for model_name, model in models:
        print(f"Fitting {model_name}")
        mf = ModelFit(model, ev_correct_by_num_ev)
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

    stmt_lkl_values = []
    labels = []
    for model_name, mf, sampler in results:
        labels.append(model_name)
        stmt_lkl_values.append(mf.stmt_err(sampler))
    plt.figure()
    plt.bar(range(len(stmt_lkl_values)), stmt_lkl_values, tick_label=labels)
    plt.ylim(bottom=250)

    # Print table of results before plotting
    table = Texttable()
    table_data = [('Model', '-log(Max Lkl) (Stmt)')]
    table_data.extend(zip(labels, stmt_lkl_values))
    table.add_rows(table_data)
    print(table.draw())

