import sys
import itertools
from matplotlib import pyplot as plt
from os.path import abspath, dirname, join
import numpy as np
from indra.tools import assemble_corpus as ac
from indra.belief import load_default_probs, SimpleScorer, BeliefEngine
from bioexp.util import pkldump, pklload
from bioexp.curation.process_curations import \
            get_correctness_data, load_curated_pkl_files, get_full_curations


def set_fitted_belief(reader_input, stmts):
    fitted_probs = {'syst': {}, 'rand': {}}
    for reader, rd_dict in reader_input.items():
        bm_file = rd_dict['belief_model']
        model_fit, sampler = pklload(bm_file)
        map_params = model_fit.get_map_params(sampler)
        fitted_probs['syst'][reader] = map_params['Syst']
        fitted_probs['rand'][reader] = map_params['Rand']
    fitted_scorer = SimpleScorer(fitted_probs)
    #fitted_scorer = MaxScorer(fitted_probs)
    be = BeliefEngine(scorer=fitted_scorer)
    be.set_prior_probs(stmts)
    be.set_hierarchy_probs(stmts)


def MaxScorer(SimpleScorer):
    def score_evidence_list(self, evidences):
        """Return belief score given a list of supporting evidences."""
        def _score(evidences):
            if not evidences:
                return 0
            # Collect all unique sources
            sources = [ev.source_api for ev in evidences]
            uniq_sources = numpy.unique(sources)
            # Calculate the systematic error factors given unique sources
            syst_factors = {s: self.prior_probs['syst'][s]
                            for s in uniq_sources}
            # Calculate the radom error factors for each source
            rand_factors = {k: [] for k in uniq_sources}
            for ev in evidences:
                rand_factors[ev.source_api].append(
                    evidence_random_noise_prior(
                        ev,
                        self.prior_probs['rand'],
                        self.subtype_probs))
            # The probability of incorrectness is the product of the
            # source-specific probabilities
            neg_prob_prior = 1
            # ########## THIS SECTION IS DIFFERENT ####################
            # Take the minimum source-specific error estimate and don't
            # multiply with other source-specific error rates
            for s in uniq_sources:
                source_specific_probs.append(
                        syst_factors[s] + numpy.prod(rand_factors[s]))
            neg_prob_prior = np.max(source_specific_probs)
            #############################################
            # Finally, the probability of correctness is one minus incorrect
            prob_prior = 1 - neg_prob_prior
            return prob_prior
        pos_evidence = [ev for ev in evidences if
                        not ev.epistemics.get('negated')]
        neg_evidence = [ev for ev in evidences if
                        ev.epistemics.get('negated')]
        pp = _score(pos_evidence)
        np = _score(neg_evidence)
        # The basic assumption is that the positive and negative evidence
        # can't simultaneously be correct.
        # There are two cases to consider. (1) If the positive evidence is
        # incorrect then there is no Statement and the belief should be 0,
        # irrespective of the negative evidence.
        # (2) If the positive evidence is correct and the negative evidence
        # is incorrect.
        # This amounts to the following formula:
        # 0 * (1-pp) + 1 * (pp * (1-np)) which we simplify below
        score = pp * (1 - np)
        return score


def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n


def plot_calibration_curve(results):
    stmts = []
    bel_tuples = []
    for bin_ix, bin_dict in results['bins'].items():
        for stmt in bin_dict['stmts']:
            stmts.append(stmt)
            correct = 0 if stmt.get_hash() in results['incorr_hashes'] else 1
            bel_tuples.append((stmt.belief, correct))
    bel_tuples.sort(key=lambda x: x[0])
    # Plot with moving average, raw
    beliefs, corrects = zip(*bel_tuples)
    beliefs = np.array(beliefs)
    corrects = np.array(corrects)
    n = 30
    plt.plot(beliefs, corrects, linestyle='', marker='.')
    plt.xlabel('Belief Score')
    plt.ylabel('Empirical Correctness')
    plt.plot(beliefs[:-(n-1)], moving_average(corrects, n=n), linestyle='',
             marker='.', color='red')
    plt.plot(beliefs, beliefs, color='gray')
    bel_pre = []
    corr_pre = []
    for bel_val, corr_group in itertools.groupby(bel_tuples,
                                                 key=lambda x: x[0]):
        corr_vals = [t[1] for t in corr_group]
        bel_pre.append(bel_val)
        corr_pre.append(np.mean(corr_vals))
    plt.plot(bel_pre[:-(n-1)], moving_average(corr_pre, n=n), marker='.',
             color='orange')

if __name__ == '__main__':
    # Prevent issues in pickling the results
    sys.setrecursionlimit(50000)

    reader_input = {
       'reach': {
         'pkl_list': [
            'bioexp_reach_sample_uncurated_19-12-14.pkl',
            'bioexp_reach_sample_uncurated_20-02-19.pkl',
            'bioexp_reach_sample_tsv.pkl'],
         'source_list': ['bioexp_paper_reach', 'bioexp_paper_tsv'],
         'belief_model':
              'reach_orig_belief_stmt_evidence_sampler'},
       'rlimsp': {
         'pkl_list': ['bioexp_rlimsp_sample_uncurated.pkl'],
         'source_list': ['bioexp_paper_rlimsp'],
         'belief_model':
              'rlimsp_orig_belief_stmt_evidence_sampler'},
       'trips': {
         'pkl_list': ['bioexp_trips_sample_uncurated.pkl'],
         'source_list': ['bioexp_paper_trips'],
         'belief_model':
              'trips_orig_belief_stmt_evidence_sampler'},
    }

    asmb_pkl = join(dirname(abspath(__file__)), '..', '..', 'data',
                    'bioexp_asmb_preassembled.pkl')
    all_stmts = ac.load_statements(asmb_pkl)
    all_stmts_by_hash = {stmt.get_hash(): stmt for stmt in all_stmts}

    curations = {}
    for reader, rd_dict in reader_input.items():
        stmts = load_curated_pkl_files(rd_dict['pkl_list'])
        stmts_dict = {stmt.get_hash(): stmt for stmt in stmts}
        curations[reader] = get_full_curations(rd_dict['source_list'],
                                            stmts_dict, aggregation='evidence')

    for r1, r2 in itertools.combinations(curations.keys(), 2):
        r1_cur = set(curations[r1].keys())
        r2_cur = set(curations[r2].keys())
        overlap_hashes = r1_cur.intersection(r2_cur)
        print(r1, r2, ": overlap ", len(overlap_hashes), " hashes")

    # For every curation hash across the 3 pickles, set aside the ones where
    # at least 1 extraction is correct; for the ones that are incorrect, check
    # if there is evidence from other readers; if so, add this to a pickle
    # to be curated

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
                                        filter_hashes=incorr_hashes_multi_src)
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

    set_fitted_belief(reader_input, curated_stmts)

    # Group curated stmts into bins
    bins = [0., 0.6, 0.8, 0.9, 0.95, 0.99, 1.0]
    stmts_by_belief = {'corr_hashes': corr_hashes,
                       'incorr_hashes': incorr_hashes,
                       'bins': {}}
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
        if n_total == 0:
            continue
        pct_corr = n_corr / n_total
        stmts_by_belief['bins'][bin_ix] = {
                'lb': lb, 'ub': ub, 'stmts': bin_stmts,
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


