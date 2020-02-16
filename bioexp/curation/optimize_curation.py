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


def add_proposed_data(correct_by_num_ev, exp_by_num_ev):
    proposed_correct_by_num_ev = deepcopy(correct_by_num_ev)
    for num_ev, nexp in exp_by_num_ev.items():
        if num_ev in correct_by_num_ev:
            p = np.mean(correct_by_num_ev[num_ev])
        else:
            p = 0.5
            proposed_correct_by_num_ev[num_ev] = []
        nexpcorr = int(np.round(p * nexp))
        proposed_correct_by_num_ev[num_ev] += \
            (nexpcorr*[1] + (nexp-nexpcorr)*[0])
    return proposed_correct_by_num_ev


if __name__ == '__main__':
    stmts = load_reach_curated_stmts()
    source_list = ('bioexp_paper_tsv', 'bioexp_paper_reach')
    correct_by_num_ev = preprocess_data(source_list, stmts)
    proposed_data = [{1: 50, 10: 50},
                     {10: 100},
                     {1: 100},
                     {1: 25, 3: 25, 6: 25, 10: 25}]
    correct_by_num_ev = {}
    for exp_by_num_ev in proposed_data:
        proposed_correct_by_num_ev = add_proposed_data(correct_by_num_ev,
                                                       exp_by_num_ev)
        samples = get_posterior_samples(proposed_correct_by_num_ev)
        pred_unc = pred_uncertainty(belief, samples, range(1, 11))
        print('Predicted overall uncertainty with %s: %.2E' %
              (str(exp_by_num_ev), pred_unc))
