import sys
from os.path import abspath, dirname, join
import json
import pickle
from collections import Counter
from bioexp.util import prefixed_file

stmts_file = join(dirname(abspath(__file__)), '..', '..', 'data',
                  'bioexp_asmb_preassembled.pkl')


with open(stmts_file, 'rb') as fh:
    print(f'Loading {stmts_file}')
    stmts = pickle.load(fh)


def get_reader_ev_pmid_distro(reader):
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

    dd = {(i+1): pmid_distro[i] for i in range(10)}
    s = sum(dd.values())
    pmid_distro_norm = {k: v/s for k, v in dd.items()}

    dd = {(i+1): ev_distro[i] for i in range(10)}
    s = sum(dd.values())
    ev_distro_norm = {k: v/s for k, v in dd.items()}

    return ev_distro_norm, pmid_distro_norm


def dump_jsons(reader):
    print(f'Dumping distributions for {reader}')
    ev_distro_norm, pmid_distro_norm = get_reader_ev_pmid_distro(reader)

    ev_file = prefixed_file(f'{reader}_stmt_evidence_distribution', 'json')
    print(f'Dumping into {ev_file}')
    with open(ev_file, 'w') as fh:
        json.dump(ev_distro_norm, fh, indent=1)

    pmid_file = prefixed_file(f'{reader}_stmt_pmid_distribution', 'json')
    print(f'Dumping into {pmid_file}')
    with open(pmid_file, 'w') as fh:
        json.dump(pmid_distro_norm, fh, indent=1)


if __name__ == '__main__':
    readers = sys.argv[1:]
    for reader in readers:
        dump_jsons(reader)
