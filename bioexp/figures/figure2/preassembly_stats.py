from bioexp.util import set_fig_params, fontsize, format_axis
import pickle
from collections import Counter
from os.path import dirname, join
from matplotlib import pyplot as plt


data_dir = join(dirname(__file__), '..', '..', '..', 'data')
build_dir = join(dirname(__file__), '..', '..', '..', 'build')

stmts_file = join(data_dir, 'bioexp_preassembled.pkl')

# Load the pickle
print("Loading statements from %s" % stmts_file)
with open(stmts_file, 'rb') as f:
    stmts = pickle.load(f)

print("%d stmts" % len(stmts))
# Get/plot evidence distribution
ev_counts = [len(s.evidence) for s in stmts]
ev_ctr = Counter(ev_counts)
ev_ctr = sorted([(k, v) for k, v in ev_ctr.items()],
                  key=lambda x: x[1], reverse=True)
"""
with open('ev_ctr.pkl', 'wb') as f:
    pickle.dump(ev_ctr, f)
with open('ev_ctr.pkl', 'rb') as f:
    ev_ctr = pickle.load(f)
"""

counts, stmts_per_count = zip(*ev_ctr)

ev_dist_fig = join(build_dir, 'fig2_evidence_distribution.pdf')

set_fig_params()
fig = plt.figure(figsize=(2, 2), dpi=150)
plt.plot(counts, stmts_per_count, linestyle='', marker='.', markersize='1')
ax = fig.gca()
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel('Mentions per statement')
ax.set_ylabel('Number of statements')
format_axis(ax, tick_padding=2)
plt.subplots_adjust(left=0.17, bottom=0.14, right=0.94, top=0.95)
plt.savefig(ev_dist_fig)
