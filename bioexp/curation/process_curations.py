import csv
import numpy
import matplotlib.pyplot as plt


def belief(num_ev, pr, ps):
    return 1 - (ps + pr ** num_ev)


if __name__ == '__main__':
    # Get a dict of all curations by UUID
    curations = {}
    with open('rlimsp_sample_curated.tsv', 'r') as fh:
        reader = csv.reader(fh, delimiter='\t')
        next(reader)
        for row in reader:
            uuid = row[2]
            correct = row[16]
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

    # Finally, calculate the mean of correctness by number of evidence
    num_evs = sorted(correct_by_num_ev.keys())
    means = [numpy.mean(correct_by_num_ev[n]) for n in num_evs]
    beliefs = [belief(n, 0.4, 0.05) for n in num_evs]
    plt.plot(num_evs, means, 'bo-', label='Empirical mean correctness')
    plt.plot(num_evs, beliefs, 'ro-', label='INDRA belief score')
    plt.ylim(0, 1)
    plt.grid(True)
    plt.xticks(num_evs)
    plt.xlabel('Number of evidence per INDRA Statement')
    plt.legend()
    plt.show()
