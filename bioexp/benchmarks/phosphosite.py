from collections import defaultdict
from bioexp.util import pklload
from indra.statements import Phosphorylation
from indra.tools import assemble_corpus as ac


def filter_stmts(stmts):
    stmts = ac.filter_by_type(stmts, Phosphorylation)
    stmts = ac.filter_genes_only(stmts, specific_only=True)
    stmts = ac.filter_human_only(stmts)
    return stmts


def stratify_by_detail(stmts):
    stmts_by_detail = defaultdict(list)
    for stmt in stmts:
        sk = stmt.sub.get_grounding()
        if stmt.enz:
            ek = stmt.enz.get_grounding()
            key = (ek, sk)
            stmts_by_detail[key].append(stmt)
            if stmt.residue:
                key = (ek, sk, stmt.residue)
                stmts_by_detail[key].append(stmt)
                if stmt.position:
                    key = (ek, sk, stmt.residue, stmt.position)
                    stmts_by_detail[key].append(stmt)
    return stmts_by_detail


def compare_by_detail(indra_stmts, psp_stmts):
    in_psp = {}
    not_in_psp = {}
    for key, stmts in indra_stmts.items():
        if key in psp_stmts:
            in_psp[key] = (stmts, psp_stmts[key])
        else:
            not_in_psp[key] = stmts
    return in_psp, not_in_psp


if __name__ == '__main__':
    indra_stmts = pklload('reading_only_asmb_preassembled')
    psp_stmts = pklload('phosphosite')

    indra_stmts = filter_stmts(indra_stmts)
    psp_stmts = filter_stmts(psp_stmts)

    indra_by_detail = stratify_by_detail(indra_stmts)
    psp_by_detail = stratify_by_detail(psp_stmts)

    in_psp, not_in_psp = compare_by_detail(indra_by_detail, psp_by_detail)