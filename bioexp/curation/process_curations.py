import sys
import csv
import numpy
import scipy.optimize
import matplotlib.pyplot as plt


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


if __name__ == '__main__':
    # Get a dict of all curations by UUID
    curation_file = sys.argv[1]

    curations = {}
    with open(curation_file, 'r') as fh:
        reader = csv.reader(fh, delimiter='\t')
        next(reader)
        for row in reader:
            uuid = row[1]
            correct = row[15]
            correct = None if correct == '' else int(correct)
            if uuid in curations:
                curations[uuid].append(correct)
            else:
                curations[uuid] = [correct]

    # Filter to only curations where every entry for the
    # given UUID was curated
    full_curations = {k: v for k, v in curations.items()
                      if all([vv is not None for vv in v])}


    # Now organize the curations by number of evidence
    correct_by_num_ev = {}
    for uuid, corrects in full_curations.items():
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
    plt.errorbar(num_evs, means, yerr=std, fmt='bo-',
                 label='Empirical mean correctness')
    plt.plot(num_evs, beliefs, 'ro-', label='INDRA belief score')
    plt.ylim(0, 1)
    plt.grid(True)
    plt.xticks(num_evs)
    plt.xlabel('Number of evidence per INDRA Statement')
    plt.legend()
    plt.show()
