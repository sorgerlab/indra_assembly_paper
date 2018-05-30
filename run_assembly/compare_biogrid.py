from indra.tools import assemble_corpus as ac
from indra.statements import Complex
from util import pklload

stmts = pklload('reading_only_asmb_preassembled')
cplx = ac.filter_by_type(stmts, Complex)
hgnc_stmts = []
for stmt in cplx:
    if all([('HGNC' in ag.db_refs) for ag in stmt.agent_list()]):
        hgnc_stmts.append(stmt)

bg_stmts = pklload('biogrid')

print("hgnc_stmts", len(hgnc_stmts))

in_biogrid = []
bg_keys = [s.matches_key() for s in bg_stmts]
for s in hgnc_stmts:
    if s.matches_key in bg_keys:
        in_biogrid.append(s)


