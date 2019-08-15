import sys
import numpy
import pickle
from copy import copy
from collections import defaultdict, Counter
from indra.tools import assemble_corpus as ac
from indra.assemblers.tsv import TsvAssembler


def sample_stmts(stmts, source, n, ev_min=1, ev_max=10):
    filt_stmts = [s for s in stmts if s.agent_list()[0] is not None]
    uuids = []
    stmts_by_uuid = defaultdict(list)
    for stmt in filt_stmts:
        for ev in stmt.evidence:
            new_stmt_flat = copy(stmt)
            if ev.source_api == source:
                uuids.append(stmt.uuid)
                new_stmt_flat.evidence = [ev]
                # Transfer the annotations from the evidence to top-level stmt
                for ix, ag in enumerate(new_stmt_flat.agent_list()):
                    if ag is None:
                        continue
                    ag.db_refs = ev.annotations['agents']['raw_grounding'][ix]
                stmts_by_uuid[stmt.uuid].append(new_stmt_flat)

    uuid_ctr = Counter(uuids)
    sampled_stmts = []
    for ev_num in range(ev_min, ev_max+1):
        filt_uuids = [uuid for uuid, count in uuid_ctr.items()
                      if count == ev_num]
        # Sample from the uuids filtered to the evidence count
        sampled_uuids = list(numpy.random.choice(filt_uuids, size=n,
                             replace=True))
        for uuid in sampled_uuids:
            for uuid_stmt in stmts_by_uuid[uuid]:
                sampled_stmts.append(uuid_stmt)
    return sampled_stmts


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

    # Sample stmts for each source
    for source in sources:
        stmts_sample = sample_stmts(stmts, source, n, ev_min, ev_max)
        make_tsv(stmts_sample, source, output_dir)
