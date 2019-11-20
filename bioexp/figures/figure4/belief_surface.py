from os.path import dirname, join
import numpy as np
from indra.statements import *
from indra.belief import BeliefEngine, SimpleScorer
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from indra.util import plot_formatting as pf
from bioexp.util import based

build_dir = based
plot_filename = join(build_dir, 'fig4_belief_surface.pdf')

ag1 = Agent('Agent1', db_refs={'HGNC':'1', 'UP':'1'})
ag2 = Agent('Agent2', db_refs={'HGNC':'2', 'UP':'2'})

example_probs = {'rand': {'source1': 0.4, 'source2': 0.3},
                 'syst': {'source1': 0.2, 'source2': 0.1}}
ss = SimpleScorer(example_probs)

be = BeliefEngine(ss)
num_ev = 5
ev1_coords = np.zeros((num_ev+1, num_ev+1))
ev2_coords = np.zeros((num_ev+1, num_ev+1))
beliefs = np.zeros((num_ev+1, num_ev+1))
for ev1_ix in range(0, num_ev+1):
    for ev2_ix in range(0, num_ev+1):
        ev1_list = [Evidence(pmid=j, source_api='source1')
                    for j in range(0, ev1_ix)]
        ev2_list = [Evidence(pmid=j, source_api='source2')
                    for j in range(0, ev2_ix+1)]
        combined_evidence = ev1_list + ev2_list
        if len(combined_evidence) == 0:
            belief=  0
        else:
            st = Phosphorylation(ag1, ag2, evidence=combined_evidence)
            be.set_prior_probs([st])
            belief = st.belief
        beliefs[ev1_ix, ev2_ix] = belief
        ev1_coords[ev1_ix, ev2_ix] = ev1_ix
        ev2_coords[ev1_ix, ev2_ix] = ev2_ix

plt.ion()
#pf.set_fig_params()
fig = plt.figure(figsize=(4, 4), dpi=150)
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(ev2_coords, ev1_coords, beliefs, cmap='plasma')
plt.show()
ax.set_xlabel('Source 2')
ax.set_ylabel('Source 1')
ax.set_zlabel('Belief')
ax.set_zlim(0.5, 1)
ax.view_init(elev=16., azim=-116.)
plt.savefig(plot_filename)
