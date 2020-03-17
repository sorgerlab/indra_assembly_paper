import json
import pickle
from collections import Counter


with open('../../data/bioexp_asmb_preassembled.pkl', 'rb') as fh:
    stmts = pickle.load(fh)


pmid_cnt = []
ev_cnt = []
for stmt in stmts:
    reach_ev = [ev for ev in stmt.evidence if ev.source_api == 'reach']
    if 1 <= len(reach_ev) <= 10:
        npmids = len({ev.pmid for ev in reach_ev})
        pmid_cnt.append(npmids) 
        ev_cnt.append(len(reach_ev))

pmid_distro = [Counter(pmid_cnt)[i] for i in range(1,11)]
ev_distro = [Counter(ev_cnt)[i] for i in range(1,11)]


with open('reach_stmt_pmid_distribution.json', 'w') as fh:
    dd = {(i+1): pmid_distro[i] for i in range(10)}
    s = sum(dd.values())
    dd = {k: v/s for k, v in dd.items()}
    json.dump(dd, fh, indent=1)

with open('stmt_evidence_distribution.json', 'w') as fh:
    dd = {(i+1): ev_distro[i] for i in range(10)}
    s = sum(dd.values())
    dd = {k: v/s for k, v in dd.items()}
    json.dump(dd, fh, indent=1)
