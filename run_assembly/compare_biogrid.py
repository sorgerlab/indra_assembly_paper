import numpy as np
from indra.tools import assemble_corpus as ac
from indra.statements import Complex
from util import pklload
from matplotlib import pyplot as plt

stmts = pklload('reading_only_preassembled')
cplx = ac.filter_by_type(stmts, Complex)
hgnc_stmts = []
for stmt in cplx:
    if all([('HGNC' in ag.db_refs) for ag in stmt.agent_list()]):
        hgnc_stmts.append(stmt)

bg_stmts = pklload('biogrid')

print("hgnc_stmts", len(hgnc_stmts))

in_biogrid = []
print("matches key 1")
bg_keys = [s.matches_key() for s in bg_stmts]
print("matches key 2")
stmt_keys = [s.matches_key() for s in hgnc_stmts]

bins = [0, 0.5, 0.8, 0.9, 0.99, 1.0]

# Sample statements from each belief bin
stmt_bins = ((0.0, 0.5, '0_0.5'),
             (0.5, 0.8, '0.5_0.8'),
             (0.8, 0.9, '0.8_0.9'),
             (0.9, 0.99, '0.9_0.99'),
             (0.99, 1.0, '0.99_1.0'))

in_bg = np.zeros(len(bins) - 1)
not_in_bg = np.zeros(len(bins) - 1)
not_in_bg_stmts = {}
for i, (lbound, ubound, label) in enumerate(stmt_bins):
    stmts_by_belief = [s for s in hgnc_stmts
                       if s.belief >= lbound and s.belief < ubound]
    all_keys = set([s.matches_key() for s in stmts_by_belief])
    in_bg_keys = all_keys.intersection(bg_keys)
    not_in_bg_keys = all_keys.difference(in_bg_keys)
    in_bg[i] = len(in_bg_keys)
    not_in_bg[i] = len(stmts_by_belief) - len(in_bg_keys)
    not_in_bg_stmts[label] = not_in_bg_keys

stmts_by_key = {s.matches_key(): s for s in hgnc_stmts}

# Plot statement counts
index = np.array(range(len(bins)-1))
plt.figure()
plt.bar(index, not_in_bg, bottom=in_bg, label='Not In BioGrid')
plt.bar(index, in_bg, label='In BioGrid')
ax = plt.gca()
ax.set_xticks(index)
ax.set_xticklabels(('< 0.5', '0.5-0.8', '0.8-0.9', '0.9-0.99',
                    '> 0.99'))
plt.ylabel('Number of Statements')
plt.xlabel('Belief scores')
plt.legend(loc='upper right')
plt.savefig('output/biogrid_stmt_numbers.pdf')

# Plot statement percentages
plt.figure()
plt.bar(index, in_bg / (in_bg + not_in_bg))
ax = plt.gca()
ax.set_xticks(index)
ax.set_xticklabels(('< 0.5', '0.5-0.8', '0.8-0.9', '0.9-0.99',
                    '> 0.99'))
plt.ylabel('Percent of Statements in BioGrid')
plt.xlabel('Belief scores')
plt.savefig('output/biogrid_stmt_percentages.pdf')


"""
for s in hgnc_stmts:
    if s.matches_key in bg_keys:
        in_biogrid.append(s)
"""

