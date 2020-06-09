import sys
from os.path import abspath, dirname, join
import json
import pickle
from collections import Counter
from bioexp.util import prefixed_file

reader = sys.argv[1]


stmts_file = join(dirname(abspath(__file__)), '..', '..', 'data',
                  'bioexp_asmb_preassembled.pkl')


with open(stmts_file, 'rb') as fh:
    stmts = pickle.load(fh)


pmid_cnt = []
ev_cnt = []
for stmt in stmts:
    reader_ev = [ev for ev in stmt.evidence if ev.source_api == reader]
    if 1 <= len(reader_ev) <= 10:
        npmids = len({ev.pmid for ev in reader_ev})
        pmid_cnt.append(npmids) 
        ev_cnt.append(len(reader_ev))


pmid_distro = [Counter(pmid_cnt)[i] for i in range(1,11)]
ev_distro = [Counter(ev_cnt)[i] for i in range(1,11)]


pmid_file = prefixed_file(f'{reader}_stmt_pmid_distribution', 'json')
with open(pmid_file, 'w') as fh:
    dd = {(i+1): pmid_distro[i] for i in range(10)}
    s = sum(dd.values())
    dd = {k: v/s for k, v in dd.items()}
    json.dump(dd, fh, indent=1)


ev_file = prefixed_file(f'{reader}_stmt_evidence_distribution', 'json')
with open(ev_file, 'w') as fh:
    dd = {(i+1): ev_distro[i] for i in range(10)}
    s = sum(dd.values())
    dd = {k: v/s for k, v in dd.items()}
    json.dump(dd, fh, indent=1)
