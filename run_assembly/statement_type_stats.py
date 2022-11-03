"""This script produces a table of Statement type statistics per reader
before normalization or filtering, i.e., at the level of raw Statements."""

import csv
from collections import defaultdict, Counter
import pickle

with open('bioexp_all_raw.pkl', 'rb') as fh:
    stmts = pickle.load(fh)

stmt_type_by_source = defaultdict(list)
for stmt in stmts:
    for ev in stmt.evidence:
        stmt_type_by_source[ev.source_api].append(stmt.__class__.__name__)

readers = ['reach', 'sparser', 'medscan', 'trips', 'rlimsp', 'isi']
all_stmt_types = set()
for reader in readers:
    all_stmt_types |= set(stmt_type_by_source[reader])

header = ['Statement type'] + readers
rows = [header]
for stmt_type in sorted(all_stmt_types):
    cnts = []
    for reader in readers:
        cnt = Counter(stmt_type_by_source[reader]).get(stmt_type, 0)
        cnts.append(str(cnt))
    rows.append([stmt_type] + cnts)

with open('reader_stmt_types.csv', 'w') as fh:
    w = csv.writer(fh)
    w.writerows(rows)
