import pickle
from util import *
from indra.mechlinker import MechLinker
from copy import deepcopy

#stmts = pklload('preassembled')
#for s in stmts:
#    s.evidence = []
#pkldump(stmts, 'preassembled_no_ev')
stmts = pklload('preassembled_no_ev')
#stmts = [s for s in stmts
#         if 'PDK1' in [a.name for a in s.agent_list() if a is not None]]
print("Deepcopying stmts")
ml_stmts = deepcopy(stmts)
ml = MechLinker(ml_stmts)

print("Gathering activities")
ml.gather_explicit_activities()
print("Reducing activities")
ml.reduce_activities()

stmts_uuid = {}
# Index statements by UUID
for s in stmts:
    stmts_uuid[s.uuid] = s

ml_stmts_uuid = {}
for s in ml.statements:
    ml_stmts_uuid[s.uuid] = s

updated = []
for uuid in stmts_uuid.keys():
    raw = stmts_uuid[uuid]
    upd = ml_stmts_uuid[uuid]
    if not raw.matches(upd):
        updated.append((raw, upd))
# 34,495 statements updated

print("Inferring mods")
inferred_mods = ml.infer_modifications(stmts)

print("Inferring regs")
inferred_regs = ml.infer_activations(stmts)

print("Inferring AFs")
inferred_afs = ml.infer_active_forms(stmts)

# Infer active forms
# Validate inferences?
