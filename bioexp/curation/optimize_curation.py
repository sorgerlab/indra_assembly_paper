import json
import logging
import numpy as np
from copy import deepcopy
from bioexp.curation.process_curations import *
from bioexp.curation.belief_models import OrigBeliefStmt
from bioexp.curation.model_fit import ModelFit, ens_sample

logger = logging.getLogger('optimize_curations')


def pred_uncertainty(fun, param_samples, x_values, x_probs=None):
    if x_probs is None:
        x_probs = np.ones((len(x_values), )) / len(x_values)
    preds = [fun(p, x_values) for p in param_samples]
    sum_var = 0
    for idx, xv in enumerate(x_values):
        xp = x_probs.get(xv, 0)
        var = np.var([p[idx] for p in preds])
        logger.info('Variance at %d evidences: %.2E' % (xv, var))
        sum_var += xp * var
    logger.info('Overall uncertainty: %.2E' % sum_var)
    return sum_var


def add_proposed_data(correct_by_num_ev, exp_by_num_ev, p):
    proposed_correct_by_num_ev = deepcopy(correct_by_num_ev)
    pexps = model.stmt_predictions(p, exp_by_num_ev)
    for idx, (num_ev, nexp) in enumerate(exp_by_num_ev.items()):
        if num_ev not in correct_by_num_ev:
            proposed_correct_by_num_ev[num_ev] = []
        nexpcorr = int(np.round(pexps[idx] * nexp))
        proposed_correct_by_num_ev[num_ev] += \
            (nexpcorr*[1] + (nexp-nexpcorr)*[0])
    return proposed_correct_by_num_ev


def curation_cost(num_ev, type='linear'):
    """Return an estimate of the cost of curation for a statement with given
    number of evidences."""
    if type == 'linear':
        return num_ev
    elif type == 'log':
        return 1 + np.log(num_ev)
    elif type == 'log2':
        return 1 + np.log2(num_ev)


def cur_for_cost(cost, num_ev, cost_type):
    """Return the number of curations that can be done at a given cost limit."""
    return int(np.floor(cost/curation_cost(num_ev, cost_type)))


def find_next_best(cost=10, maxev=10, ev_probs=None):
    # First, establish reference values based on current curations
    ref_sampler = ens_sample(mf,
                             nwalkers=default_nwalkers,
                             burn_steps=1000,
                             sample_steps=default_nsteps,
                             )
    ref_samples = ref_sampler.flatchain
    ref_pred_unc = pred_uncertainty(mf.model.stmt_predictions,
                                    ref_samples,
                                    range(1, maxev+1),
                                    x_probs=ev_probs)
    logger.info('Baseline prediction uncertainty: %.2E' % ref_pred_unc)
    ref_param_unc = np.var(ref_samples, 0)
    logger.info('Baseline parameter variance: er: %.2E es: %.2E' %
                tuple(ref_param_unc))
    ref_params = np.mean(ref_samples, 0)
    logger.info('Baseline parameter means: er: %.2E es: %.2E' %
                tuple(ref_params))

    deltas = {}

    for num_ev in range(1, maxev+1):
        # Propose experiments with num_ev evidences
        proposed_curations_by_num_ev = {num_ev: cur_for_cost(cost, num_ev,
                                                             default_cost_type)}
        proposed_correct_by_num_ev = \
            add_proposed_data(mf.data,
                              proposed_curations_by_num_ev,
                              ref_params)
        mf_prop = ModelFit(model, proposed_correct_by_num_ev)
        prop_sampler = ens_sample(mf_prop,
                                  nwalkers=default_nwalkers,
                                  burn_steps=1000,
                                  sample_steps=default_nsteps,
                                  )
        samples = prop_sampler.flatchain
        pred_unc = pred_uncertainty(mf_prop.model.stmt_predictions,
                                    samples,
                                    range(1, maxev+1),
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
    logger.info('Delta summary:')
    for num_ev, delta in deltas_sorted:
        logger.info('%s: %.2E' % (num_ev, delta))
    opt_num_ev = deltas_sorted[0][0]
    logger.info('You should next curate %d statements with %d evidences.' %
                (cur_for_cost(cost, opt_num_ev, default_cost_type),
                 opt_num_ev))
    return opt_num_ev, cur_for_cost(cost, opt_num_ev, default_cost_type)


if __name__ == '__main__':
    default_nsteps = 10000
    default_nwalkers = 50
    default_cost_type = 'log2'
    stmts = load_reader_curated_stmts(reader='rlimsp')
    ev_probs = load_stmt_evidence_distribution()
    source_list = ['bioexp_paper_rlimsp']
    model = OrigBeliefStmt()
    stmt_correct_by_num_ev = get_correctness_data(source_list, stmts,
                                                  aggregation='evidence')
    mf = ModelFit(model, stmt_correct_by_num_ev)
    opt_i, num_i = find_next_best(cost=50, maxev=10, ev_probs=ev_probs)
