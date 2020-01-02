import sys
import numpy
import pickle
from copy import copy
from collections import defaultdict, Counter
from indra.tools import assemble_corpus as ac
from indra.assemblers.tsv import TsvAssembler
from bioexp.util import pkldump



def sample_stmts(stmts, source, n, ev_min=1, ev_max=10):
    filt_stmts = [s for s in stmts if s.agent_list()[0] is not None]
    uuids = []
    stmts_by_uuid_tsv = defaultdict(list)
    stmts_by_uuid_pkl = {}
    # Create new synthetic statements containing only evidence from the
    # specified source
    for stmt in filt_stmts:
        new_ev_list = []
        for ev in stmt.evidence:
            new_stmt_flat = copy(stmt)
            if ev.source_api == source:
                new_ev_list.append(ev)
                uuids.append(stmt.uuid)
                new_stmt_flat.evidence = [ev]
                # Transfer the annotations from the evidence to top-level stmt
                for ix, ag in enumerate(new_stmt_flat.agent_list()):
                    if ag is None:
                        continue
                    ag.db_refs = ev.annotations['agents']['raw_grounding'][ix]
                    ag.db_refs['TEXT'] = \
                        ev.annotations['agents']['raw_text'][ix]
                stmts_by_uuid_tsv[stmt.uuid].append(new_stmt_flat)
        if new_ev_list:
            new_stmt = copy(stmt)
            new_stmt.evidence = new_ev_list
            stmts_by_uuid_pkl[stmt.uuid] = new_stmt

    assert set(uuids) == set(stmts_by_uuid_tsv.keys())
    assert set(uuids) == set(stmts_by_uuid_pkl.keys())
    uuid_ctr = Counter(uuids)
    sampled_stmts_tsv = []
    sampled_stmts_pkl = []

    for ev_num in range(ev_min, ev_max+1):
        filt_uuids = [uuid for uuid, count in uuid_ctr.items()
                      if count == ev_num]
        if not filt_uuids:
            continue
        # Sample from the uuids filtered to the evidence count
        # If the number of uuids with this evidence count is less than or
        # equal to the number for sampling, simply take all of them
        if len(filt_uuids) <= n:
            sampled_uuids = filt_uuids
        # Otherwise, sample with replacement
        else:
            sampled_uuids = list(numpy.random.choice(filt_uuids, size=n,
                                 replace=True))
        for uuid in sampled_uuids:
            # Add statements with multiple evidences to list for pkl
            sampled_stmts_pkl.append(stmts_by_uuid_pkl[uuid])
            # Add all flattened (1 evidence) statements to list for tsv
            for uuid_stmt in stmts_by_uuid_tsv[uuid]:
                sampled_stmts_tsv.append(uuid_stmt)
    sampled_uuids = [s.uuid for s in sampled_stmts_pkl]
    print(len(sampled_uuids), len(set(sampled_uuids)))
    return sampled_stmts_tsv, sampled_stmts_pkl


def make_tsv(stmts, source, output_dir):
    ta = TsvAssembler(stmts)
    output_file = '%s/%s_sample_uncurated.tsv' % (output_dir, source)
    res = ta.make_model(output_file, add_curation_cols=True)


if __name__ == '__main__':
    stmts_pkl = sys.argv[1]
    n = int(sys.argv[2])
    ev_min = int(sys.argv[3])
    ev_max = int(sys.argv[4])
    output_dir = sys.argv[5]
    sources = sys.argv[6:]

    stmts = ac.load_statements(stmts_pkl)

    # Set numpy random seed
    numpy.random.seed(1)

    sys.setrecursionlimit(50000)

    # Sample stmts for each source
    for source in sources:
        print(source)
        stmts_sample_tsv, stmts_sample_pkl = sample_stmts(stmts, source, n,
                                                          ev_min, ev_max)
        make_tsv(stmts_sample_tsv, source, output_dir)
        pkldump(stmts_sample_pkl, '%s_sample_uncurated' % source)

