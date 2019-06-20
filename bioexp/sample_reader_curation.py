import sys
import numpy
import pickle
from indra.assemblers.tsv import TsvAssembler


def sample_stmts(stmts, source, n):
    counts = [len([ev for ev in stmt.evidence if ev.source_api == source])
              for stmt in stmts]
    sum_counts = sum(counts)
    probs = [c/sum_counts for c in counts]
    sampled_stmts = list(numpy.random.choice(stmts, size=n, replace=True,
                                             p=probs))
    return sampled_stmts


def make_tsv(stmts, source):
    ta = TsvAssembler(stmts)
    output_file = '%s_sample.tsv' % source
    res = ta.make_model(output_file, add_curation_cols=True)


if __name__ == '__main__':
    stmts_pkl = sys.argv[1]
    source = sys.argv[2]
    n = int(sys.argv[3])

    with open(stmts_pkl, 'rb') as fh:
        stmts = pickle.load(fh)

    stmts_sample = sample_stmts(stmts, source, n)
    make_tsv(stmts_sample, source)