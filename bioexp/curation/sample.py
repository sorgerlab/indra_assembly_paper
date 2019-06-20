import sys
import numpy
import pickle
from copy import copy
from indra.tools import assemble_corpus as ac
from indra.assemblers.tsv import TsvAssembler

def sample_stmts(stmts, source, n):
    asmb_stmts = []
    flat_stmts = []
    for stmt in stmts:
        new_stmt_asmb = copy(stmt)
        new_stmt_flat = copy(stmt)
        new_ev_list = []
        for ev in stmt.evidence:
            if ev.source_api == source:
                new_ev_list.append(ev)
                new_stmt_flat.evidence = [ev]
                # Transfer the annotations from the evidence to top-level stmt
                for ix, ag in enumerate(new_stmt_flat.agent_list()):
                    if ag is None:
                        continue
                    ag.db_refs = ev.annotations['agents']['raw_grounding'][ix]
                flat_stmts.append(new_stmt_flat)
        if new_ev_list:
            new_stmt_asmb.evidence = new_ev_list
            asmb_stmts.append(new_stmt_asmb)

    """
    counts = [len([ev for ev in stmt.evidence for stmt in source_stmts]
    sum_counts = sum(counts)
    probs = [c/sum_counts for c in counts]
    #sampled_stmts = list(numpy.random.choice(stmts, size=n, replace=True,
    #                                         p=probs))
    """
    sampled_stmts = list(numpy.random.choice(flat_stmts, size=n, replace=True))
    return sampled_stmts


def make_tsv(stmts, source, output_dir):
    ta = TsvAssembler(stmts)
    output_file = '%s/%s_sample_uncurated.tsv' % (output_dir, source)
    res = ta.make_model(output_file, add_curation_cols=True)


if __name__ == '__main__':
    stmts_pkl = sys.argv[1]
    n = int(sys.argv[2])
    output_dir = sys.argv[3]
    sources = sys.argv[4:]

    stmts = ac.load_statements(stmts_pkl)

    # Sample stmts for each source
    for source in sources:
        stmts_sample = sample_stmts(stmts, source, n)
        make_tsv(stmts_sample, source, output_dir)
