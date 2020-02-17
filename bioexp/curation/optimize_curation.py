from copy import deepcopy
from bioexp.curation.process_curations import *


def pred_uncertainty(fun, param_samples, x_values, x_probs=None):
    if x_probs is None:
        x_probs = np.ones((len(x_values), )) / len(x_values)
    sum_var = 0
    for xv, xp in zip(x_values, x_probs):
        var = xp * np.var([fun(xv, *p) for p in param_samples])
        print('%d: %.2E' % (xv, var))
        sum_var += var
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


def proposal_uncertainty(exp_by_num_ev):
    proposed_correct_by_num_ev = add_proposed_data(correct_by_num_ev,
                                                   exp_by_num_ev,
                                                   (0.3, 0.3))
    samples = get_posterior_samples(proposed_correct_by_num_ev)
    pred_unc = pred_uncertainty(belief, samples, range(1, 11))
    print('Predicted overall uncertainty with %s: %.2E' %
          (str(exp_by_num_ev), pred_unc))
    return pred_unc


def optimize_proposal_uncertainty_simple(cost=100):
    # NOTE: due to the discrete nature of the optimization problem,
    # this mode of optimization, which relies on Jacobian evaluations
    # does not work well in practice.
    bounds = [(0, cost)] * 10
    cost_constraint = {'type': 'ineq',
                       'fun': lambda x: (sum(x * range(1, len(x)+1)) - cost)}
    fun = lambda x: proposal_uncertainty({i+1: int(np.floor(x[i]))
                                          for i in range(10)})
    res = scipy.optimize.minimize(fun,
                                  bounds=bounds,
                                  constraints=[cost_constraint],
                                  x0=([cost] + [0]*9))
    return res.x


def optimize_proposal_uncertainty_anneal(cost=100):
    # Both parameters have to be between 0 and 1
    cost_constraint = lambda x: all(sum(x * range(1, len(x)+1)) < cost)
    positive_constraint = lambda x: all(x >= 0)
    accept_test = lambda x: positive_constraint(x) and cost_constraint(x)
    fun = lambda x: proposal_uncertainty({i+1: int(np.floor(x[i]))
                                          for i in range(10)})
    res = scipy.optimize.basinhopping(fun,
                                      accept_test=accept_test,
                                      x0=([cost] + [0]*9))
    return res.x


def find_next_best(cost=10):
    samples = get_posterior_samples(correct_by_num_ev)
    ref_pred_unc = pred_uncertainty(belief, samples, range(1, 11))
    deltas = {}
    for i in range(1, 11):
        exp_by_num_ev = {i: int(np.floor(cost/i))}
        proposed_correct_by_num_ev = \
            add_proposed_data(correct_by_num_ev, exp_by_num_ev,
                              (0.3, 0.3))
        samples = get_posterior_samples(proposed_correct_by_num_ev)
        pred_unc = pred_uncertainty(belief, samples, range(1, 11))
        delta_pred_unc = ref_pred_unc - pred_unc
        deltas[i] = delta_pred_unc
        print('Curating %d statements with %d evidences '
              'will decrease belief uncertainty by %.2E.' %
              (exp_by_num_ev[i], i, delta_pred_unc))
    opt_i = sorted(deltas.items(), key=lambda x: x[1],
                   reverse=True)[0][0]
    return opt_i, int(np.floor(cost/opt_i))


if __name__ == '__main__':
    stmts = load_reach_curated_stmts()
    source_list = ('bioexp_paper_tsv', 'bioexp_paper_reach')
    correct_by_num_ev = preprocess_data(source_list, stmts)
    opt_i, num_i = find_next_best(20)
