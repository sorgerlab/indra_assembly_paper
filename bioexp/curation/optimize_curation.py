import json
import logging
from copy import deepcopy
from bioexp.curation.process_curations import *

logger = logging.getLogger('optimize_curations')


def pred_uncertainty(fun, param_samples, x_values, x_probs=None):
    if x_probs is None:
        x_probs = np.ones((len(x_values), )) / len(x_values)
    sum_var = 0
    for xv, xp in zip(x_values, x_probs):
        var = np.var([fun(xv, *p) for p in param_samples])
        logger.info('Variance at %d evidences: %.2E' % (xv, var))
        sum_var += xp * var
    return sum_var


def add_proposed_data(correct_by_num_ev, exp_by_num_ev, p):
    proposed_correct_by_num_ev = deepcopy(correct_by_num_ev)
    for num_ev, nexp in exp_by_num_ev.items():
        if num_ev not in correct_by_num_ev:
            proposed_correct_by_num_ev[num_ev] = []
        pexp = belief(num_ev, *p)
        nexpcorr = int(np.round(pexp * nexp))
        proposed_correct_by_num_ev[num_ev] += \
            (nexpcorr*[1] + (nexp-nexpcorr)*[0])
    return proposed_correct_by_num_ev


def proposal_uncertainty(exp_by_num_ev, assumed_params,
                         maxev=10, ev_probs=None):
    proposed_correct_by_num_ev = add_proposed_data(correct_by_num_ev,
                                                   exp_by_num_ev,
                                                   assumed_params)
    samples = get_posterior_samples(proposed_correct_by_num_ev,
                                    nsteps=10000)
    pred_unc = pred_uncertainty(belief, samples, range(1, maxev+1),
                                x_probs=ev_probs)
    logger.info('Predicted overall uncertainty with %s: %.2E' %
                (str(exp_by_num_ev), pred_unc))
    return pred_unc


def optimize_proposal_uncertainty_simple(cost=100, maxev=10):
    # NOTE: due to the discrete nature of the optimization problem,
    # this mode of optimization, which relies on Jacobian evaluations
    # does not work well in practice.

    # Establish a baseline based on current data
    samples = get_posterior_samples(correct_by_num_ev,
                                    nsteps=10000)
    ref_params = np.mean(samples, 0)

    bounds = [(0, cost)] * maxev
    cost_constraint = {'type': 'ineq',
                       'fun': lambda x: (sum(x * range(1, len(x)+1)) - cost)}
    fun = lambda x: proposal_uncertainty({i+1: int(np.floor(x[i]))
                                          for i in range(maxev)},
                                         ref_params, maxev=10,
                                         ev_probs=ev_probs)
    res = scipy.optimize.minimize(fun,
                                  bounds=bounds,
                                  constraints=[cost_constraint],
                                  x0=([cost] + [0]*(maxev-1)))
    return res.x


def optimize_proposal_uncertainty_anneal(cost=100, maxev=10):
    # Establish a baseline based on current data
    samples = get_posterior_samples(correct_by_num_ev,
                                    nsteps=10000)
    ref_params = np.mean(samples, 0)

    # Both parameters have to be between 0 and 1
    cost_constraint = lambda x: sum(x * range(1, len(x)+1)) < cost
    positive_constraint = lambda x: all(x >= 0)
    accept_test = lambda x: positive_constraint(x) and cost_constraint(x)
    fun = lambda x: proposal_uncertainty({i+1: int(np.floor(x[i]))
                                          for i in range(maxev)},
                                         ref_params, maxev=10,
                                         ev_probs=ev_probs)
    res = scipy.optimize.basinhopping(fun,
                                      accept_test=accept_test,
                                      x0=([cost] + [0]*maxev))
    return res.x


def find_next_best(cost=10, maxev=10, ev_probs=None):
    # First, establish reference values based on current curations
    samples = get_posterior_samples(correct_by_num_ev, nsteps=10000)
    ref_pred_unc = pred_uncertainty(belief, samples, range(1, maxev+1),
                                    x_probs=ev_probs)
    ref_param_unc = np.var(samples, 0)
    ref_params = np.mean(samples, 0)

    deltas = {}

    def cur_for_num_ev(cost, num_ev):
        """Return the number of curations that can be done at a given cost"""
        return int(np.floor(cost/num_ev))

    for num_ev in range(1, maxev+1):
        # Propose experiments with i evidences
        proposed_curations_by_num_ev = {num_ev: cur_for_num_ev(cost, num_ev)}
        proposed_correct_by_num_ev = \
            add_proposed_data(correct_by_num_ev, proposed_curations_by_num_ev,
                              ref_params)
        samples = get_posterior_samples(proposed_correct_by_num_ev,
                                        nsteps=10000)
        pred_unc = pred_uncertainty(belief, samples, range(1, maxev+1),
                                    x_probs=ev_probs)
        logger.info('Predicted overall uncertainty with %s: %.2E' %
                    (str(proposed_curations_by_num_ev), pred_unc))

        # Calculate decrease in prediction uncertainty
        delta_pred_unc = ref_pred_unc - pred_unc
        deltas[num_ev] = delta_pred_unc
        logger.info('Curating %d statements with %d evidences '
                    'will decrease belief uncertainty by %.2E.' %
                    (proposed_curations_by_num_ev[num_ev], num_ev,
                     delta_pred_unc))
        # Calculate delta in parameter uncertainty
        param_unc = np.var(samples, 0)
        delta_param_unc = ref_param_unc - param_unc
        logger.info('It will also decrease parameter uncertainty by '
                    'er: %.2E and es: %.2E' % tuple(delta_param_unc))
    deltas_sorted = sorted(deltas.items(), key=lambda x: x[1],
                           reverse=True)
    logger.info('Deltas per num evidence: %s' % str(deltas_sorted))
    opt_num_ev = deltas_sorted[0][0]
    logger.info('You should next curate %d statements with %d evidences.' %
                (cur_for_num_ev(cost, opt_num_ev), opt_num_ev))
    return opt_num_ev, cur_for_num_ev(cost, opt_num_ev)


if __name__ == '__main__':
    with open('stmt_evidence_distribution.json', 'r') as fh:
        ev_probs = json.load(fh)
        ev_probs = {int(k): v for k, v in ev_probs.items()}
    stmts = load_reach_curated_stmts()
    source_list = ('bioexp_paper_tsv', 'bioexp_paper_reach')
    correct_by_num_ev = preprocess_data(source_list, stmts)
    opt_i, num_i = find_next_best(20, ev_probs=ev_probs)
