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


def optimize_proposal_uncertainty(cost=100):
    # Both parameters have to be between 0 and 1
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


if __name__ == '__main__':
    stmts = load_reach_curated_stmts()
    source_list = ('bioexp_paper_tsv', 'bioexp_paper_reach')
    correct_by_num_ev = preprocess_data(source_list, stmts)
    optimize_proposal_uncertainty()
